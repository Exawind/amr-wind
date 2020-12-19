#include "amr-wind/equation_systems/vof/SplitAdvection.H"
#include "amr-wind/equation_systems/vof/lagrangian_advection.H"
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
    amrex::Real dt)
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
    p += fluxR.size();

    if (isweep % 3 == 0) {
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real wL = wmac(i, j, k);
                amrex::Real wR = wmac(i, j, k + 1);
                multiphase::lagrangian_advection(
                    i, j, k, 3, dtdz, wL, wR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 3, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
                    domhi.z);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real uL = umac(i, j, k);
                amrex::Real uR = umac(i + 1, j, k);
                multiphase::lagrangian_advection(
                    i, j, k, 1, dtdx, uL, uR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 1, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
                    domhi.x);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real vL = vmac(i, j, k);
                amrex::Real vR = vmac(i, j + 1, k);
                multiphase::lagrangian_advection(
                    i, j, k, 2, dtdy, vL, vR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 2, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
                    domhi.y);
            });
    } else if (isweep % 2 == 0) {
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real vL = vmac(i, j, k);
                amrex::Real vR = vmac(i, j + 1, k);
                multiphase::lagrangian_advection(
                    i, j, k, 2, dtdy, vL, vR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 2, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
                    domhi.y);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real wL = wmac(i, j, k);
                amrex::Real wR = wmac(i, j, k + 1);
                multiphase::lagrangian_advection(
                    i, j, k, 3, dtdz, wL, wR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 3, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
                    domhi.z);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real uL = umac(i, j, k);
                amrex::Real uR = umac(i + 1, j, k);
                multiphase::lagrangian_advection(
                    i, j, k, 1, dtdx, uL, uR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 1, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
                    domhi.x);
            });
    } else {
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real uL = umac(i, j, k);
                amrex::Real uR = umac(i + 1, j, k);
                multiphase::lagrangian_advection(
                    i, j, k, 1, dtdx, uL, uR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 1, volfrac, fluxL, fluxC, fluxR, pbc, domlo.x,
                    domhi.x);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real vL = vmac(i, j, k);
                amrex::Real vR = vmac(i, j + 1, k);
                multiphase::lagrangian_advection(
                    i, j, k, 2, dtdy, vL, vR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 2, volfrac, fluxL, fluxC, fluxR, pbc, domlo.y,
                    domhi.y);
            });
        amrex::ParallelFor(
            bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real wL = wmac(i, j, k);
                amrex::Real wR = wmac(i, j, k + 1);
                multiphase::lagrangian_advection(
                    i, j, k, 3, dtdz, wL, wR, volfrac, fluxL, fluxC, fluxR);
            });
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                multiphase::balance_fluxes(
                    i, j, k, 3, volfrac, fluxL, fluxC, fluxR, pbc, domlo.z,
                    domhi.z);
            });
    }

    // Remove VOF debris
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        multiphase::remove_vof_debris(i, j, k, volfrac);
    });
}

} // namespace amr_wind