#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLCodedInlet.H"
#include "amr-wind/wind_energy/ABLFillCodedInlet.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include <AMReX_PlotFileUtil.H>

#include <sstream>
#include <iostream>
#include <string>
#include <dlfcn.h>

namespace amr_wind {

ABLCodedInlet::ABLCodedInlet(CFDSim& sim)
    : m_sim(sim)
    , m_time(m_sim.time())
    , m_repo(m_sim.repo())
    , m_mesh(m_sim.mesh())
    , m_velocity(m_sim.repo().get_field("velocity"))
    , m_temperature(m_sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("codedInlet");

    pp.query("lib", m_user_lib);
    if (m_user_lib != "") {
        if (FILE *lib = fopen(m_user_lib.c_str(), "r")) {
            // lib found
            fclose(lib);

            // symbols from shared library
            userfun_lib = dlopen(m_user_lib.c_str(), RTLD_NOW);
            m_user_coded.get_vel = reinterpret_cast<Vfun_handle>(
                    dlsym(userfun_lib, "velocity"));
            m_user_coded.get_temp = reinterpret_cast<Tfun_handle>(
                    dlsym(userfun_lib, "temperature"));

            // attempt to use user-coded functions
            double Vtmp[3];
            m_user_coded.get_vel(0,0,0,0,Vtmp);
            amrex::Print() << "ABLCodedInlet: Loaded xvelocity_field function "
                << "V(0,0,0,0)=" << Vtmp[0] << " " << Vtmp[1] << " " << Vtmp[2] << std::endl;
            double Ttmp;
            m_user_coded.get_temp(0,0,0,0,Ttmp);
            amrex::Print() << "ABLCodedInlet: Loaded temperature_field function "
                << "T(0,0,0,0)=" << Ttmp << std::endl;

            m_active = true;

        } else {
            amrex::Print() << "ABLCodedInlet: Shared library not found "
                << m_user_lib << std::endl;
        }
    }
}

ABLCodedInlet::~ABLCodedInlet()
{
    if (m_active) {
        dlclose(userfun_lib);
    }
}

void ABLCodedInlet::post_init_actions()
{
    if (m_active) {
        m_velocity.register_fill_patch_op<ABLFillCodedInlet>(m_mesh, m_time, *this);
        m_temperature.register_fill_patch_op<ABLFillCodedInlet>(m_mesh, m_time, *this);
    }
}

void ABLCodedInlet::pre_advance_work() {}

void ABLCodedInlet::post_advance_work() {}

void ABLCodedInlet::set_velocity(
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

    BL_PROFILE("amr-wind::ABLCodedInlet::set_velocity");

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
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

            amrex::ParallelFor(
                bx, [=,coded=m_user_coded] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    double Vtmp[3];
                    coded.get_vel(time,x,y,z,Vtmp);

                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vels = {
                        AMREX_D_DECL(Vtmp[0],Vtmp[1],Vtmp[2])};
                    for (int n = 0; n < numcomp; n++) {
                        arr(i, j, k, dcomp + n) = vels[orig_comp + n];
                    }
                });
        }
    }
}

void ABLCodedInlet::set_temperature(
    const int lev,
    const amrex::Real time,
    const Field& fld,
    amrex::MultiFab& mfab) const
{
    if (!m_active) {
        return;
    }

    BL_PROFILE("amr-wind::ABLCodedInlet::set_temperature");

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

            amrex::ParallelFor(
                bx, [=,coded=m_user_coded] AMREX_GPU_DEVICE(
                        int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                    coded.get_temp(time,x,y,z,arr(i,j,k));
                });
        }
    }
}

} // namespace amr_wind