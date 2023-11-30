#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLUserDefined.H"
#include "amr-wind/wind_energy/ABLFillUDF.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include <AMReX_PlotFileUtil.H>

#include <sstream>
#include <iostream>
#include <string>
#include <dlfcn.h>

namespace amr_wind {

ABLUserDefined::ABLUserDefined(CFDSim& sim)
    : m_sim(sim)
    , m_time(m_sim.time())
    , m_repo(m_sim.repo())
    , m_mesh(m_sim.mesh())
    , m_velocity(m_sim.repo().get_field("velocity"))
    , m_temperature(m_sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("UDF");

#ifdef AMR_WIND_USE_EXTERNAL_UDF
    // load external UDF library
    pp.query("inflow_lib", m_inflow_lib);
    if (!m_inflow_lib.empty()) {
        if (FILE* lib = fopen(m_inflow_lib.c_str(), "r")) {
            // lib found
            fclose(lib);

            // symbols from shared library
            userfun_lib = dlopen(m_inflow_lib.c_str(), RTLD_NOW);
            m_udf.get_vel =
                reinterpret_cast<Vfun_ptr>(dlsym(userfun_lib, "velocity"));
            m_udf.get_temp =
                reinterpret_cast<Tfun_ptr>(dlsym(userfun_lib, "temperature"));

            // RTLD_NOW does not seem to guarantee that symbols are loaded
            // immediately -- do a test function call here to verify the user
            // code is well behaved.
            double Vtmp[3];
            m_udf.get_vel(0, 0, 0, 0, Vtmp);
            amrex::Print() << "ABLUserDefined: Loaded xvelocity_field function"
                           << std::endl;
            double Ttmp;
            m_udf.get_temp(0, 0, 0, 0, Ttmp);
            amrex::Print()
                << "ABLUserDefined: Loaded temperature_field function"
                << std::endl;

            m_active = true;

        } else {
            amrex::Print() << "ABLUserDefined: Shared library not found "
                           << m_inflow_lib << std::endl;
        }
    }
#endif
}

ABLUserDefined::~ABLUserDefined()
{
    if (m_active) {
        dlclose(userfun_lib);
    }
}

void ABLUserDefined::post_init_actions()
{
    if (m_active) {
        m_velocity.register_fill_patch_op<ABLFillUDF>(
            m_mesh, m_time, *this);
        m_temperature.register_fill_patch_op<ABLFillUDF>(
            m_mesh, m_time, *this);
    }
}

void ABLUserDefined::pre_advance_work() {}

void ABLUserDefined::post_advance_work() {}

void ABLUserDefined::set_velocity(
    const int lev,
    const amrex::Real time,
    const Field& fld,
    amrex::MultiFab& mfab,
    const int dcomp,
    const int orig_comp) const
{
    if (!m_active) {
        return;
    }

    BL_PROFILE("amr-wind::ABLUserDefined::set_velocity");

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (bctype[ori] != BC::mass_inflow) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            if (!gbx.cellCentered()) {
                gbx.enclosedCells();
            }
            const auto& bx = gbx & dbx;
            if (!bx.ok()) {
                continue;
            }

            const auto& arr = mfab[mfi].array();
            const int numcomp = mfab.nComp();

#ifdef AMR_WIND_USE_EXTERNAL_UDF
            amrex::ParallelFor(
                bx, [=, udf = m_udf] AMREX_GPU_DEVICE(
                        int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    double Vtmp[3];
                    udf.get_vel(time, x, y, z, Vtmp);

                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vels = {
                        AMREX_D_DECL(Vtmp[0], Vtmp[1], Vtmp[2])};
                    for (int n = 0; n < numcomp; n++) {
                        arr(i, j, k, dcomp + n) = vels[orig_comp + n];
                    }
                });
#endif
        }
    }
}

void ABLUserDefined::set_temperature(
    const int lev,
    const amrex::Real time,
    const Field& fld,
    amrex::MultiFab& mfab) const
{
    if (!m_active) {
        return;
    }

    BL_PROFILE("amr-wind::ABLUserDefined::set_temperature");

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (bctype[ori] != BC::mass_inflow) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            if (!gbx.cellCentered()) {
                gbx.enclosedCells();
            }
            const auto& bx = gbx & dbx;
            if (!bx.ok()) {
                continue;
            }

            const auto& arr = mfab[mfi].array();

#ifdef AMR_WIND_USE_EXTERNAL_UDF
            amrex::ParallelFor(
                bx, [=, udf = m_udf] AMREX_GPU_DEVICE(
                        int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    udf.get_temp(time, x, y, z, arr(i, j, k));
                });
#endif
        }
    }
}

} // namespace amr_wind
