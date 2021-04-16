#include "amr-wind/incflo.H"
#include "amr-wind/diffusion/diffusion.H"

using namespace amrex;

namespace diffusion {

Vector<Array<LinOpBCType, AMREX_SPACEDIM>> get_diffuse_tensor_bc(
    amr_wind::Field& velocity, Orientation::Side side) noexcept
{
    const auto& geom = velocity.repo().mesh().Geom(0);
    Vector<Array<LinOpBCType, AMREX_SPACEDIM>> r(3);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom.isPeriodic(dir)) {
            r[0][dir] = LinOpBCType::Periodic;
            r[1][dir] = LinOpBCType::Periodic;
            r[2][dir] = LinOpBCType::Periodic;
        } else {
            auto bc = velocity.bc_type()[Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::zero_gradient: {
                // All three components are Neumann
                r[0][dir] = LinOpBCType::Neumann;
                r[1][dir] = LinOpBCType::Neumann;
                r[2][dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall: {
                // All three components are Dirichlet
                r[0][dir] = LinOpBCType::Dirichlet;
                r[1][dir] = LinOpBCType::Dirichlet;
                r[2][dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::slip_wall: {
                // Tangential components are Neumann
                // Normal     component  is  Dirichlet
                r[0][dir] = LinOpBCType::Neumann;
                r[1][dir] = LinOpBCType::Neumann;
                r[2][dir] = LinOpBCType::Neumann;

                r[dir][dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::wall_model: {
                // Tangential components are inhomogeneous Neumann
                // Normal     component  is  Dirichlet
                r[0][dir] = LinOpBCType::inhomogNeumann;
                r[1][dir] = LinOpBCType::inhomogNeumann;
                r[2][dir] = LinOpBCType::inhomogNeumann;

                r[dir][dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}

Array<LinOpBCType, AMREX_SPACEDIM>
get_diffuse_scalar_bc(amr_wind::Field& scalar, Orientation::Side side) noexcept
{
    Array<LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (scalar.repo().mesh().Geom(0).isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = scalar.bc_type()[Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::zero_gradient:
            case BC::slip_wall: {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::wall_model:
            case BC::fixed_gradient: {
                r[dir] = LinOpBCType::inhomogNeumann;
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall: {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}

Array<MultiFab, AMREX_SPACEDIM> average_velocity_eta_to_faces(
    const amrex::Geometry& geom, MultiFab const& cc_eta)
{
    BL_PROFILE("amr-wind::diffusion::average_velocity_eta_to_faces");
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    Array<MultiFab, AMREX_SPACEDIM> r{
        {MultiFab(
             amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1, 0,
             MFInfo(), fact),
         MultiFab(
             amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1, 0,
             MFInfo(), fact),
         MultiFab(
             amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1, 0,
             MFInfo(), fact)}};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, geom);
    fixup_eta_on_domain_faces(geom, r, cc_eta);
    return r;
}

void fixup_eta_on_domain_faces(
    const amrex::Geometry& geom,
    Array<MultiFab, AMREX_SPACEDIM>& fc,
    MultiFab const& cc)
{
    BL_PROFILE("amr-wind::diffusion::fixup_eta_on_domain_faces");

    const Box& domain = geom.Domain();
    MFItInfo mfi_info{};
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cc, mfi_info); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const& cca = cc.const_array(mfi);

        int idim = 0;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i - 1, j, k);
                    });
            }
        }

        idim = 1;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i, j - 1, k);
                    });
            }
        }

        idim = 2;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) = cca(i, j, k - 1);
                    });
            }
        }
    }
}
} // namespace diffusion
