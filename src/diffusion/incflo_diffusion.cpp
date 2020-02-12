#include <incflo.H>

using namespace amrex;

bool extrapolate=false;// fixme make an input maybe? extrapolate could cause negative viscosity if not used carefully

DiffusionTensorOp*
incflo::get_diffusion_tensor_op ()
{
    if (!m_diffusion_tensor_op) m_diffusion_tensor_op.reset(new DiffusionTensorOp(this));
    return m_diffusion_tensor_op.get();
}

DiffusionScalarOp*
incflo::get_diffusion_scalar_op ()
{
    if (!m_diffusion_scalar_op) m_diffusion_scalar_op.reset(new DiffusionScalarOp(this));
    return m_diffusion_scalar_op.get();
}

Vector<Array<LinOpBCType,AMREX_SPACEDIM> >
incflo::get_diffuse_tensor_bc (Orientation::Side side) const noexcept
{
    Vector<Array<LinOpBCType,AMREX_SPACEDIM>> r(3);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (Geom(0).isPeriodic(dir)) {
            r[0][dir] = LinOpBCType::Periodic;
            r[1][dir] = LinOpBCType::Periodic;
            r[2][dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
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

Array<LinOpBCType,AMREX_SPACEDIM>
incflo::get_diffuse_scalar_bc (Orientation::Side side) const noexcept
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (Geom(0).isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = m_bc_type[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            case BC::slip_wall:
            case BC::no_slip_wall:
            case BC::wall_model: //fixme this should be an inhomogNeumann for wall_model and slip_wall for heat flux
            {
                r[dir] = LinOpBCType::Neumann;
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

void
incflo::wall_model_bc(int lev, amrex::Real utau, amrex::Real umag, const amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM>& fc_eta, amrex::MultiFab& density, amrex::MultiFab& velocity) const
{

    const Geometry& gm = Geom(lev);
    const Box& domain = gm.Domain();
    MFItInfo mfi_info{};

    // copies cell center to face
    Real c0 = 1.0;
    Real c1 = 0.0;
    
    // linear extrapolate onto face
    if(extrapolate)
    {
        c0 =  1.5;
        c1 = -0.5;
    }
    
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(density,mfi_info); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto den = density.array(mfi);

        int idim = 0;
//        if (!gm.isPeriodic(idim)) {
//            Array4<Real> const& fca = fc[idim].array(mfi);
//            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
//                amrex::ParallelFor(amrex::bdryLo(bx, idim),
//                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i+1,j,k);
//                });
//            }
//            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
//                amrex::ParallelFor(amrex::bdryHi(bx, idim),
//                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                    fca(i,j,k) = c0*cca(i-1,j,k) + c1*cca(i-2,j,k);
//                });
//            }
//        }
//
//        idim = 1;
//        if (!gm.isPeriodic(idim)) {
//            Array4<Real> const& fca = fc[idim].array(mfi);
//            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
//                amrex::ParallelFor(amrex::bdryLo(bx, idim),
//                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i,j+1,k);
//
//                });
//            }
//            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
//                amrex::ParallelFor(amrex::bdryHi(bx, idim),
//                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                    fca(i,j,k) = c0*cca(i,j-1,k) + c1*cca(i,j-2,k);
//
//                });
//            }
//        }

        idim = 2;
        if (!gm.isPeriodic(idim)) {
            Array4<Real const> const& eta = fc_eta[idim]->const_array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // density and velocity are cell centered
                    // viscosity eta is face centered
                    Real rho = c0*den(i,j,k) + c1*den(i,j,k+1);
                    Real mu = eta(i,j,k);
                    Real vx = c0*vel(i,j,k,0)+ c1*vel(i,j,k+1,0);
                    Real vy = c0*vel(i,j,k,1)+ c1*vel(i,j,k+1,1);
                    
                    // for convience velocity also holds BC's which is why derivatives are going into the ghost cells at k-1
                    // fixme this is confusing, maybe have a separate mfab for BCs?
                    
                    // dudz(x,y,z=0)
                    vel(i,j,k-1,0) = rho*utau*utau*vx/umag/mu;
                    // dvdz(x,y,z=0)
                    vel(i,j,k-1,1) = rho*utau*utau*vy/umag/mu;
                    // w(x,y,z=0)
                    vel(i,j,k-1,2) = 0.0;
                    
                    //fixme remove this print
                    //if(i==0 and j==0) amrex::Print() << "wall model at 0,0. dudz = " << vel(i,j,k-1,0) << " dvdz= " << vel(i,j,k-1,1) << " vx/umag= " << vx/umag << " vy/umag= "  << vy/umag << " mu= "  << mu << " utau= " << utau <<  std::endl;

                });
            }
//            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
//                amrex::ParallelFor(amrex::bdryHi(bx, idim),
//                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//                {
//                    fca(i,j,k) = c0*cca(i,j,k-1) + c1*cca(i,j,k-2);
//
//                });
//            }
        }
    }
    
}

Array<MultiFab,AMREX_SPACEDIM>
incflo::average_velocity_eta_to_faces (int lev, MultiFab const& cc_eta) const
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
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, Geom(lev));
    fixup_eta_on_domain_faces(lev, r, cc_eta);
    return r;
}

Array<MultiFab,AMREX_SPACEDIM>
incflo::average_tracer_eta_to_faces (int lev, int comp, MultiFab const& cc_eta) const
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
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc, Geom(lev));
    fixup_eta_on_domain_faces(lev, r, cc);
    return r;
}

void
incflo::fixup_eta_on_domain_faces (int lev, Array<MultiFab,AMREX_SPACEDIM>& fc,
                                   MultiFab const& cc) const
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
    
    const Geometry& gm = Geom(lev);
    const Box& domain = gm.Domain();
    MFItInfo mfi_info{};
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cc,mfi_info); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const& cca = cc.const_array(mfi);

        int idim = 0;
        if (!gm.isPeriodic(idim)) {
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
        if (!gm.isPeriodic(idim)) {
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
        if (!gm.isPeriodic(idim)) {
            Array4<Real> const& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j,k) + c1*cca(i,j,k+1);
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fca(i,j,k) = c0*cca(i,j,k-1) + c1*cca(i,j,k-2);
                });
            }
        }
    }
}
