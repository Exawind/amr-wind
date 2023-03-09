#ifndef HYDROSTATIC_OPS_H
#define HYDROSTATIC_OPS_H
#include <AMReX_MultiFabUtil.H>
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::hydrostatic {

static void define_rho0(
    amr_wind::Field& rho0,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Vector<amrex::Geometry> geom)
{
    for (int lev = 0; lev < rho0.repo().num_active_levels(); ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        for (amrex::MFIter mfi(rho0(lev)); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.validbox();
            auto rho0_arr = rho0(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real zbtm = problo[2] + k * dx[2];
                    const amrex::Real vof =
                        amrex::max(amrex::min(1.0, (wlev - zbtm) / dx[2]), 0.0);
                    rho0_arr(i, j, k) = vof * rho1 + (1.0 - vof) * rho2;
                });
        }
    }
}

static void define_p0(
    amr_wind::Field& p0,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Real grav_z,
    const amrex::Vector<amrex::Geometry> geom)
{
    for (int lev = 0; lev < p0.repo().num_active_levels(); ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();
        for (amrex::MFIter mfi(p0(lev)); mfi.isValid(); ++mfi) {
            amrex::Box const& nbx = mfi.grownnodaltilebox();
            auto p0_arr = p0(lev).array(mfi);
            amrex::ParallelFor(
                nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Height of pressure node
                    const amrex::Real hnode = k * dx[2];
                    // Liquid height
                    const amrex::Real hliq = wlev - problo[2];
                    // Integrated (top-down in z) phase heights to pressure node
                    amrex::Real ih_g = amrex::max(
                        0.0, amrex::min(probhi[2] - hliq, probhi[2] - hnode));
                    amrex::Real ih_l = amrex::max(
                        0.0, amrex::min(hliq - hnode, hliq - problo[2]));
                    // Integrated rho at pressure node
                    const amrex::Real irho = rho1 * ih_l + rho2 * ih_g;

                    // Add term to reference pressure
                    p0_arr(i, j, k) = -irho * grav_z;
                });
        }
    }
}

} // namespace amr_wind::hydrostatic

#endif