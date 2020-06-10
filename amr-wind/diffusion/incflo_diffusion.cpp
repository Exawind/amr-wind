#include "amr-wind/incflo.H"
#include "amr-wind/diffusion/diffusion.H"

using namespace amrex;

bool extrapolate=false;// fixme make an input maybe? extrapolate could cause negative viscosity if not used carefully

namespace diffusion {

Vector<Array<LinOpBCType, AMREX_SPACEDIM>>
get_diffuse_tensor_bc(amr_wind::Field& velocity, Orientation::Side side) noexcept
{
    const auto& geom = velocity.repo().mesh().Geom(0);
    Vector<Array<LinOpBCType,AMREX_SPACEDIM>> r(3);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom.isPeriodic(dir)) {
            r[0][dir] = LinOpBCType::Periodic;
            r[1][dir] = LinOpBCType::Periodic;
            r[2][dir] = LinOpBCType::Periodic;
        } else {
            auto bc = velocity.bc_type()[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::zero_gradient:
            {
                // All three components are Neumann
                r[0][dir] = LinOpBCType::Neumann;
                r[1][dir] = LinOpBCType::Neumann;
                r[2][dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall:
            {
                // All three components are Dirichlet
                r[0][dir] = LinOpBCType::Dirichlet;
                r[1][dir] = LinOpBCType::Dirichlet;
                r[2][dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::slip_wall:
            {
                // Tangential components are Neumann
                // Normal     component  is  Dirichlet
                r[0][dir] = LinOpBCType::Neumann;
                r[1][dir] = LinOpBCType::Neumann;
                r[2][dir] = LinOpBCType::Neumann;

                r[dir][dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::wall_model:
            {
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
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (scalar.repo().mesh().Geom(0).isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = scalar.bc_type()[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::zero_gradient:
            case BC::slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::wall_model:
            case BC::fixed_gradient:
            {
                r[dir] = LinOpBCType::inhomogNeumann;
                break;
            }
            case BC::mass_inflow:
            case BC::no_slip_wall:
            {
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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real shear_stress(
    int i,
    int j,
    int k,
    amrex::Real utau2,
    amrex::Real umag,
    amrex::Array4<amrex::Real const> const& rho,
    amrex::Array4<amrex::Real const> const& vel,
    int comp) noexcept
{
    return rho(i, j, k) * utau2 * vel(i, j, k, comp) / umag;
}


void wall_model_bc(
    amr_wind::Field& velocity,
    const amrex::Real utau,
    const amrex::Real umag,
    const amr_wind::FieldState fstate)
{
    BL_PROFILE("amr-wind::diffusion::wall_model_bc");
    auto& repo = velocity.repo();
    auto& density = repo.get_field("density", fstate);
    auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();

    amrex::Orientation xlo(amrex::Direction::x, amrex::Orientation::low);
    amrex::Orientation ylo(amrex::Direction::y, amrex::Orientation::low);
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation xhi(amrex::Direction::x, amrex::Orientation::high);
    amrex::Orientation yhi(amrex::Direction::y, amrex::Orientation::high);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);

    // copies cell center to face
    Real c0 = 1.0;
    Real c1 = 0.0;

    // linear extrapolate onto face
    if(extrapolate)
    {
        c0 =  1.5;
        c1 = -0.5;
    }

    const Real utau2 = utau*utau;

    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        MFItInfo mfi_info{};

        auto& rho_lev = density(lev);
        auto& vel_lev = velocity(lev);
        auto& eta_lev = viscosity(lev);

        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            auto vel = vel_lev.array(mfi);
            auto bc  = vel_lev.array(mfi);
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            int idim = 0;

            if (!geom.isPeriodic(idim)) {
                if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                    velocity.bc_type()[xlo] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i, j, k) + c1 * eta(i + 1, j, k);
                            // Dirichlet BC
                            bc(i - 1, j, k, 0) = 0.0;
                            // Inhomogeneous Neumann BC
                            bc(i - 1, j, k, 1) = shear_stress(i, j, k, utau2, umag, den, vel, 1) / mu;
                            bc(i - 1, j, k, 2) = shear_stress(i, j, k, utau2, umag, den, vel, 2) / mu;
                        });
                }

                if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                    velocity.bc_type()[xhi] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryHi(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i - 1, j, k) + c1 * eta(i - 2, j, k);
                            // Dirichlet BC's
                            bc(i, j, k, 0) = 0.0;
                            // Inhomogeneous Neumann BC
                            bc(i, j, k, 1) = shear_stress(i - 1, j, k, utau2, umag, den, vel, 1) / mu;
                            bc(i, j, k, 2) = shear_stress(i - 1, j, k, utau2, umag, den, vel, 2) / mu;
                        });
                }
            }

            idim = 1;

            if (!geom.isPeriodic(idim)) {
                if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                    velocity.bc_type()[ylo] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i, j, k) + c1 * eta(i, j + 1, k);
                            // Inhomogeneous Neumann BC
                            bc(i, j - 1, k, 0) = shear_stress(i, j, k, utau2, umag, den, vel, 0) / mu;
                            // Dirichlet BC
                            bc(i, j - 1, k, 1) = 0.0;
                            // Inhomogeneous Neumann BC
                            bc(i, j - 1, k, 2) = shear_stress(i, j, k, utau2, umag, den, vel, 2) / mu;
                        });
                }

                if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                    velocity.bc_type()[yhi] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryHi(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i, j - 1, k) + c1 * eta(i, j - 2, k);
                            // Inhomogeneous Neumann BC
                            bc(i, j, k, 0) = shear_stress(i, j - 1, k, utau2, umag, den, vel, 0) / mu;
                            // Dirichlet BC
                            bc(i, j, k, 1) = 0.0;
                            // Inhomogeneous Neumann BC
                            bc(i, j, k, 2) = shear_stress(i, j - 1, k, utau2, umag, den, vel, 2) / mu;
                        });
                }
            }

            idim = 2;

            if (!geom.isPeriodic(idim)) {
                if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                    velocity.bc_type()[zlo] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i, j, k) + c1 * eta(i, j, k + 1);
                            // Inhomogeneous Neumann BC
                            bc(i, j, k - 1, 0) = shear_stress(i, j, k, utau2, umag, den, vel, 0) / mu;
                            bc(i, j, k - 1, 1) = shear_stress(i, j, k, utau2, umag, den, vel, 1) / mu;
                            // Dirichlet BC
                            bc(i, j, k - 1, 2) = 0.0;
                        });
                }

                if (bx.bigEnd(idim) == domain.bigEnd(idim) &&
                    velocity.bc_type()[zhi] ==  BC::wall_model) {
                    amrex::ParallelFor(
                        amrex::bdryHi(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const Real mu  = c0 * eta(i, j, k - 1) + c1 * eta(i, j, k - 2);
                            // Inhomogeneous Neumann BC
                            bc(i, j, k, 0) = shear_stress(i, j, k - 1, utau2, umag, den, vel, 0) / mu;
                            bc(i, j, k, 1) = shear_stress(i, j, k - 1, utau2, umag, den, vel, 1) / mu;
                            // Dirichlet BC
                            bc(i, j, k, 2) = 0.0;
                        });
                }
            }

        } // MFIter loop
    } // level loop
}

Array<MultiFab,AMREX_SPACEDIM>
average_velocity_eta_to_faces (const amrex::Geometry& geom, MultiFab const& cc_eta)
{
    BL_PROFILE("amr-wind::diffusion::average_velocity_eta_to_faces");
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    Array<MultiFab,AMREX_SPACEDIM> r{{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)}};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, geom);
    fixup_eta_on_domain_faces(geom, r, cc_eta);
    return r;
}

Array<MultiFab,AMREX_SPACEDIM>
average_tracer_eta_to_faces (int comp, const amrex::Geometry& geom, MultiFab const& cc_eta)
{
    BL_PROFILE("amr-wind::diffusion::average_tracer_eta_to_faces");
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    MultiFab cc(cc_eta, amrex::make_alias, comp, 1);
    Array<MultiFab,AMREX_SPACEDIM> r{{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)}};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc, geom);
    fixup_eta_on_domain_faces(geom, r, cc);
    return r;
}

void
fixup_eta_on_domain_faces (const amrex::Geometry& geom, Array<MultiFab,AMREX_SPACEDIM>& fc,
                                   MultiFab const& cc)
{
    BL_PROFILE("amr-wind::diffusion::fixup_eta_on_domain_faces");
    // copies cell center to face
    Real c0 = 1.0;
    Real c1 = 0.0;
    
    // linear extrapolate onto face
    if(extrapolate)
    {
        c0 =  1.5;
        c1 = -0.5;
    }
    
    const Box& domain = geom.Domain();
    MFItInfo mfi_info{};
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cc,mfi_info); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const& cca = cc.const_array(mfi);

        int idim = 0;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i+1,j,k);
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i-1,j,k) + c1*cca(i-2,j,k);
                });
            }
        }

        idim = 1;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i,j+1,k);

                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j-1,k) + c1*cca(i,j-2,k);
                });
            }
        }

        idim = 2;
        if (!geom.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i,j,k+1);
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        fca(i, j, k) =
                            c0 * cca(i, j, k - 1) + c1 * cca(i, j, k - 2);
                    });
            }
        }
    }
}
} // namespace diffusion
