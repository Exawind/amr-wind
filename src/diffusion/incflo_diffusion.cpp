#include "incflo.H"
#include "diffusion.H"

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
            case BC::no_slip_wall:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            case BC::slip_wall:
            case BC::wall_model:
            {
                r[dir] = LinOpBCType::inhomogNeumann;
                break;
            }
            case BC::mass_inflow:
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

void wall_model_bc(
    amr_wind::Field& velocity,
    const amrex::Real utau,
    const amrex::Real umag,
    const amr_wind::FieldState fstate)
{
    auto& repo = velocity.repo();
    auto& density = repo.get_field("density", fstate);
    auto& viscosity = repo.get_field("velocity_nueff");
    const int nlevels = repo.num_active_levels();

    // Wall model hard coded to be only in the zlo direction
    const int idim = 2;
    amrex::Orientation olo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation ohi(amrex::Direction::z, amrex::Orientation::high);
    AMREX_ALWAYS_ASSERT(velocity.bc_type()[olo] ==  BC::wall_model);
    AMREX_ALWAYS_ASSERT(velocity.bc_type()[ohi] != BC::wall_model);

    // copies cell center to face
    Real c0 = 1.0;
    Real c1 = 0.0;

    // linear extrapolate onto face
    if(extrapolate)
    {
        c0 =  1.5;
        c1 = -0.5;
    }

    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        AMREX_ALWAYS_ASSERT(!geom.isPeriodic(idim));
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
            auto den = rho_lev.array(mfi);
            auto eta = eta_lev.array(mfi);

            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // density and velocity are cell centered
                        // viscosity eta is face centered
                        Real rho =
                            c0 * den(i, j, k) + c1 * den(i, j, k + 1);
                        Real mu = c0 * eta(i, j, k) + c1 * eta(i, j, k + 1);
                        Real vx =
                            c0 * vel(i, j, k, 0) + c1 * vel(i, j, k + 1, 0);
                        Real vy =
                            c0 * vel(i, j, k, 1) + c1 * vel(i, j, k + 1, 1);

                        // inhomogeneous Neumann BC's
                        // mu dudz = rho utau^2
                        // dudz(x,y,z=0)
                        vel(i, j, k - 1, 0) =
                            rho * utau * utau * vx / umag / mu;
                        // dvdz(x,y,z=0)
                        vel(i, j, k - 1, 1) =
                            rho * utau * utau * vy / umag / mu;

                        // Dirichlet BC's
                        // w(x,y,z=0)
                        vel(i, j, k - 1, 2) = 0.0;
                    });
            }
        }
    }
}

void heat_flux_bc(amr_wind::Field& scalar)
{
    AMREX_ALWAYS_ASSERT(scalar.num_comp() == 1);
    const int nlevels = scalar.repo().num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        heat_flux_model_bc(lev, scalar, 0);
    }
}

void
heat_flux_model_bc(const int lev, amr_wind::Field& scalar, const int comp)
{

    const Geometry& geom = scalar.repo().mesh().Geom(lev);
    const Box& domain = geom.Domain();
    MFItInfo mfi_info{};
 
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(scalar(lev),mfi_info); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real> const& bc_a = scalar(lev).array(mfi);
        int idim = 0;

        // fixme this assume periodic
        if (!geom.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[0][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdx
                    bc_a(i-1,j,k) = local_m_bc_tracer_d;
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[3][comp];
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdx
                    bc_a(i,j,k) = local_m_bc_tracer_d ;
                });
            }
        }

        idim = 1;
        if (!geom.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[1][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdy
                    bc_a(i,j-1,k) = local_m_bc_tracer_d;
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[4][comp];
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdy
                    bc_a(i,j,k) = local_m_bc_tracer_d;
                });
            }
        }

        idim = 2;
        if (!geom.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[2][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdz
                    bc_a(i,j,k-1) = local_m_bc_tracer_d;
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = scalar.bc_values_device()[5][comp];
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdz
                    bc_a(i,j,k) = local_m_bc_tracer_d;
                });
            }
        }
    }
}


Array<MultiFab,AMREX_SPACEDIM>
average_velocity_eta_to_faces (const amrex::Geometry& geom, MultiFab const& cc_eta)
{
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    Array<MultiFab,AMREX_SPACEDIM> r{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, geom);
    fixup_eta_on_domain_faces(geom, r, cc_eta);
    return r;
}

Array<MultiFab,AMREX_SPACEDIM>
average_tracer_eta_to_faces (int comp, const amrex::Geometry& geom, MultiFab const& cc_eta)
{
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    MultiFab cc(cc_eta, amrex::make_alias, comp, 1);
    Array<MultiFab,AMREX_SPACEDIM> r{MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                              dm, 1, 0, MFInfo(), fact),
                                     MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                              dm, 1, 0, MFInfo(), fact)};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc, geom);
    fixup_eta_on_domain_faces(geom, r, cc);
    return r;
}

void
fixup_eta_on_domain_faces (const amrex::Geometry& geom, Array<MultiFab,AMREX_SPACEDIM>& fc,
                                   MultiFab const& cc)
{
    
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
