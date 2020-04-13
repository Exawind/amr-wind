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
            auto bc = velocity().bc_type()[Orientation(dir,side)];
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
            auto bc = tracer().bc_type()[Orientation(dir,side)];
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

void
incflo::wall_model_bc(const int lev,
                      const amrex::Real utau,
                      const amrex::Real umag,
                      const amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>& fc_eta,
                      const amrex::MultiFab& density,
                      const amrex::MultiFab& velocity,
                      amrex::MultiFab& bc) const
{

    amrex::Print() << "warning wall model being called with hard coded bc's, wall model is assumed on bottom and slip on top" << std::endl;

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
        const Box& bx = mfi.validbox();
        auto vel = velocity.array(mfi);
        auto den = density.array(mfi);
        auto bc_a = bc.array(mfi);
        
        int idim = 0;
        // fixme this assume periodic
        if (!gm.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int , int , int ) noexcept
                {
                    amrex::Abort("wall model bc assumes periodic should not be in here xlo");
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int , int , int ) noexcept
                {
                    amrex::Abort("wall model bc assumes periodic should not be in here xhi");
                });
            }
        }

        idim = 1;
        if (!gm.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int , int , int ) noexcept
                {
                    amrex::Abort("wall model bc assumes periodic should not be in here ylo");
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int , int , int ) noexcept
                {
                    amrex::Abort("wall model bc assumes periodic should not be in here yhi");
                });
            }
        }

            
        idim = 2;
        if (!gm.isPeriodic(idim)) {
            Array4<Real const> const& eta = fc_eta[idim]->const_array(mfi);
            
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                // fixme tried to grow box in i and j but eta is not defined on corners :(
//                amrex::Print() << amrex::bdryLo(bx, idim) << std::endl;
//                auto bxc = Box(amrex::bdryLo(bx, idim)).grow(IntVect(1,1,0));
//                amrex::Print() << bxc << std::endl;

                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {

                    // density and velocity are cell centered
                    // viscosity eta is face centered
                    Real rho = c0*den(i,j,k) + c1*den(i,j,k+1);
                    Real mu = eta(i,j,k);
                    Real vx = c0*vel(i,j,k,0)+ c1*vel(i,j,k+1,0);
                    Real vy = c0*vel(i,j,k,1)+ c1*vel(i,j,k+1,1);
          
                    // inhomogeneous Neumann BC's
                    // mu dudz = rho utau^2
                    // dudz(x,y,z=0)
                    bc_a(i,j,k-1,0) = rho*utau*utau*vx/umag/mu;
                    // dvdz(x,y,z=0)
                    bc_a(i,j,k-1,1) = rho*utau*utau*vy/umag/mu;
                    
                    // Dirichlet BC's
                    // w(x,y,z=0)
                    bc_a(i,j,k-1,2) = 0.0;
              


                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    
                    // bc_a(i,j,k,0) = 0.0; // Neumann does not need a bc since it assumes dudz=0
                    // bc_a(i,j,k,1) = 0.0; // Neumann does not need a bc since it assumes dvdz=0
                    
                    // w(x,y,z=L) = 0
                    bc_a(i,j,k,2) = 0.0; // set dirichlet on top bc
    
              
                });
            }
        }
    }
    
}


void
incflo::heat_flux_model_bc(const int lev, const int comp, amrex::MultiFab& bc) const
{

    const Geometry& gm = Geom(lev);
    const Box& domain = gm.Domain();
    MFItInfo mfi_info{};
 
    // fixme still hard coding that potential temperature is the first tracer
    AMREX_ALWAYS_ASSERT(m_ntrac==1);
    
    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(bc,mfi_info); mfi.isValid(); ++mfi) {
        
        Box const& bx = mfi.validbox();
        Array4<Real> const& bc_a = bc.array(mfi);
        
        int idim = 0;
        
        // fixme this assume periodic
        if (!gm.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[0][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdx
                    bc_a(i-1,j,k) = local_m_bc_tracer_d;
                    
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[3][comp];
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdx
                    bc_a(i,j,k) = local_m_bc_tracer_d ;
                    
                });
            }
        }

        idim = 1;
        if (!gm.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[1][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdy
                    bc_a(i,j-1,k) = local_m_bc_tracer_d;
                    
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[4][comp];
                amrex::ParallelFor(amrex::bdryHi(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdy
                    bc_a(i,j,k) = local_m_bc_tracer_d;
                    
                });
            }
        }

        idim = 2;
        if (!gm.isPeriodic(idim)) {
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[2][comp];
                amrex::ParallelFor(amrex::bdryLo(bx, idim),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // inhomogeneous Neumann BC's dTdz
                    bc_a(i,j,k-1) = local_m_bc_tracer_d;
                    
                });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                const Real local_m_bc_tracer_d = tracer().bc_values_device()[5][comp];
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
