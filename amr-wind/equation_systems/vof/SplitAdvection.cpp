
#include "amr-wind/equation_systems/vof/SplitAdvection.H"
#include "amr-wind/equation_systems/vof/split_advection.H"
#include <AMReX_Geometry.H>

using namespace amrex;

namespace amr_wind {

void multiphase::split_advection(
    int lev,
    amrex::Box const& bx,
    int isweep,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* pbc,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt,
    bool is_lagrangian)
{
    BL_PROFILE("amr-wind::multiphase::split_advection");
    Box const& bxg1 = amrex::grow(bx, 1);

    const Real dx = geom[lev].CellSize(0);
    const Real dy = geom[lev].CellSize(1);
    const Real dz = geom[lev].CellSize(2);
    Real dtdx = dt / dx;
    Real dtdy = dt / dy;
    Real dtdz = dt / dz;

    Box const& domain = geom[lev].Domain();
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    Array4<Real> fluxL = makeArray4(p, bxg1, 1);
    p += fluxL.size();
    Array4<Real> fluxC = makeArray4(p, bxg1, 1);
    p += fluxC.size();
    Array4<Real> fluxR = makeArray4(p, bxg1, 1);
    p += fluxR.size(); // NOLINT: Value not read warning

    if (!is_lagrangian) {
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::c_mask(i, j, k, volfrac, fluxC);
            });
    }

    if (isweep % 3 == 0) {
        sweep(
            2, bx, dtdz, wmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
            domhi.z, is_lagrangian);
        sweep(
            0, bx, dtdx, umac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
            domhi.x, is_lagrangian);
        sweep(
            1, bx, dtdy, vmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
            domhi.y, is_lagrangian);
    } else if (isweep % 3 == 1) {
        sweep(
            1, bx, dtdy, vmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
            domhi.y, is_lagrangian);
        sweep(
            2, bx, dtdz, wmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
            domhi.z, is_lagrangian);
        sweep(
            0, bx, dtdx, umac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
            domhi.x, is_lagrangian);
    } else {
        sweep(
            0, bx, dtdx, umac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
            domhi.x, is_lagrangian);
        sweep(
            1, bx, dtdy, vmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
            domhi.y, is_lagrangian);
        sweep(
            2, bx, dtdz, wmac, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
            domhi.z, is_lagrangian);
    }

    // Remove VOF debris
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        multiphase::remove_vof_debris(i, j, k, volfrac);
    });
}

void multiphase::sweep(
    const int dir,
    amrex::Box const& bx,
    amrex::Real dtdx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxL,
    amrex::Array4<amrex::Real> const& fluxC,
    amrex::Array4<amrex::Real> const& fluxR,
    amrex::BCRec const* pbc,
    const int dimLow,
    const int dimHigh,
    bool is_lagrangian)
{

    BL_PROFILE("amr-wind::multiphase::sweep");

    Box const& bxg1 = amrex::grow(bx, 1);

    int ii = (dir == 0) ? 1 : 0;
    int jj = (dir == 1) ? 1 : 0;
    int kk = (dir == 2) ? 1 : 0;

    if (is_lagrangian) {
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real velL = vel_mac(i, j, k);
                amrex::Real velR = vel_mac(i + ii, j + jj, k + kk);
                multiphase::lagrangian_explicit(
                    i, j, k, dir, dtdx, velL, velR, volfrac, fluxL, fluxC,
                    fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_lagrangian_fluxes(
                    i, j, k, dir, volfrac, fluxL, fluxC, fluxR, pbc, dimLow,
                    dimHigh);
            });
    } else {
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real velL = vel_mac(i, j, k);
                amrex::Real velR = vel_mac(i + ii, j + jj, k + kk);
                multiphase::eulerian_implicit(
                    i, j, k, dir, dtdx, velL, velR, volfrac, fluxL, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real velL = vel_mac(i, j, k) * dtdx;
                amrex::Real velR = vel_mac(i + ii, j + jj, k + kk) * dtdx;
                multiphase::balance_eulerian_fluxes(
                    i, j, k, dir, velL, velR, volfrac, fluxL, fluxC, fluxR, pbc,
                    dimLow, dimHigh);
            });
    }
}

} // namespace amr_wind
