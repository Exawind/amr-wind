#include "Godunov.H"
#include "MOL.H"

using namespace amrex;

void
godunov::compute_convective_term (amr_wind::FieldRepo& repo,
                                  const amr_wind::FieldState fstate,
                                  Real dt,
                                  bool constant_density,
                                  bool advect_tracer,
                                  bool godunov_ppm)
{

    // single state fields
    auto& u_mac = repo.get_field("u_mac");
    auto& v_mac = repo.get_field("v_mac");
    auto& w_mac = repo.get_field("w_mac");
    auto& vel_forces = repo.get_field("velocity_src_term");
    auto& tra_forces = repo.get_field("temperature_src_term");

    // state dependent fields
    auto& vel = repo.get_field("velocity",fstate);//fixme this is godunov so always old state?
    auto& den = repo.get_field("density",fstate);
    auto& tra = repo.get_field("temperature",fstate);
    auto& conv_u = repo.get_field("velocity_conv_term");
    auto& conv_r = repo.get_field("conv_density",fstate);
    auto& conv_t = repo.get_field("temperature_conv_term");

    const int ntrac = tra.num_comp();
    auto& geom = repo.mesh().Geom();

    // fixme moved from init_advection
    // need to add ability to change with input file
    // could also save this somewhere else instead of creating each time
    amrex::Gpu::DeviceVector<int> iconserv_velocity_d;
    amrex::Gpu::DeviceVector<int> iconserv_density_d;
    amrex::Gpu::DeviceVector<int> iconserv_tracer_d;
    iconserv_velocity_d.resize(AMREX_SPACEDIM, 0);
    iconserv_density_d.resize(1, 1);
    iconserv_tracer_d.resize(ntrac, 1);


    for (int lev = 0; lev < repo.num_active_levels(); ++lev)
    {
        u_mac(lev).FillBoundary(geom[lev].periodicity());
        v_mac(lev).FillBoundary(geom[lev].periodicity());
        w_mac(lev).FillBoundary(geom[lev].periodicity());

        MFItInfo mfi_info;
        // if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,16,16)).SetDynamic(true);
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,1024,1024)).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(den(lev),mfi_info); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();

            auto rho_arr = den(lev).array(mfi);
            auto tra_arr = tra(lev).array(mfi);
            
            // fixme would it be easier to just make a scratch rhotrac?
            Box rhotrac_box = amrex::grow(bx,2);
            rhotrac_box.grow(1);

            FArrayBox rhotracfab;
            Array4<Real> rhotrac;
            if (advect_tracer) {
                rhotracfab.resize(rhotrac_box, ntrac);
                rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho_arr(i,j,k) * tra_arr(i,j,k,n);
                });
            }

            int nmaxcomp = AMREX_SPACEDIM;
            if (advect_tracer) nmaxcomp = std::max(nmaxcomp, ntrac);

            FArrayBox tmpfab(amrex::grow(bx,1), nmaxcomp*14+1);

            godunov::compute_advection(lev, bx, AMREX_SPACEDIM,
                                      conv_u(lev).array(mfi),
                                      vel(lev).const_array(mfi),
                                      u_mac(lev).const_array(mfi),
                                      v_mac(lev).const_array(mfi),
                                      w_mac(lev).const_array(mfi),
                                      vel_forces(lev).const_array(mfi),
                                      vel.bcrec_device().data(),
                                      iconserv_velocity_d.data(),
                                      tmpfab.dataPtr(),
                                      geom, dt, godunov_ppm);

            if(!constant_density) {
                godunov::compute_advection(lev, bx, 1,
                                      conv_r(lev).array(mfi),
                                      den(lev).const_array(mfi),
                                      u_mac(lev).const_array(mfi),
                                      v_mac(lev).const_array(mfi),
                                      w_mac(lev).const_array(mfi),
                                      {},
                                      den.bcrec_device().data(),
                                      iconserv_density_d.data(),
                                      tmpfab.dataPtr(),
                                      geom, dt, godunov_ppm);
            }
            if (advect_tracer) {
                godunov::compute_advection(lev, bx, ntrac,
                                          conv_t(lev).array(mfi),
                                          rhotrac,
                                          u_mac(lev).const_array(mfi),
                                          v_mac(lev).const_array(mfi),
                                          w_mac(lev).const_array(mfi),
                                          tra_forces(lev).const_array(mfi),
                                          tra.bcrec_device().data(),
                                          iconserv_tracer_d.data(),
                                          tmpfab.dataPtr(),
                                          geom, dt, godunov_ppm);
            }

            Gpu::streamSynchronize();

        }
    }
}

void
mol::compute_convective_term (amr_wind::FieldRepo& repo,
                                 const amr_wind::FieldState fstate,
                                 bool constant_density,
                                 bool advect_tracer)
{


    // single state fields
    auto& u_mac = repo.get_field("u_mac");
    auto& v_mac = repo.get_field("v_mac");
    auto& w_mac = repo.get_field("w_mac");

    // state dependent fields
    auto& vel = repo.get_field("velocity",fstate);//fixme this is godunov so always old state?
    auto& den = repo.get_field("density",fstate);
    auto& tra = repo.get_field("temperature",fstate);
    auto& conv_u = repo.get_field("velocity_conv_term",fstate);
    auto& conv_r = repo.get_field("conv_density",fstate);
    auto& conv_t = repo.get_field("temperature_conv_term",fstate);

    const int ntrac = tra.num_comp();
    auto& geom = repo.mesh().Geom();

    int nmaxcomp = AMREX_SPACEDIM;
    if (advect_tracer) nmaxcomp = std::max(nmaxcomp,ntrac);

    for (int lev = 0; lev < repo.num_active_levels(); ++lev)
    {
        MFItInfo mfi_info;
        // if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,16,16)).SetDynamic(true);
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling(IntVect(1024,1024,1024)).SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(den(lev),mfi_info); mfi.isValid(); ++mfi)
        {

            Box const& bx = mfi.tilebox();
            auto rho_arr = den(lev).array(mfi);
            auto tra_arr = tra(lev).array(mfi);

            Box rhotrac_box = amrex::grow(bx,2);
            FArrayBox rhotracfab;
            Elixir eli_rt;
            Array4<Real> rhotrac;
            if (advect_tracer) {
                rhotracfab.resize(rhotrac_box, ntrac);
                eli_rt = rhotracfab.elixir();
                rhotrac = rhotracfab.array();
                amrex::ParallelFor(rhotrac_box, ntrac,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    rhotrac(i,j,k,n) = rho_arr(i,j,k) * tra_arr(i,j,k,n);
                });
            }

            {
                Box tmpbox = amrex::surroundingNodes(bx);
                int tmpcomp = nmaxcomp*AMREX_SPACEDIM;

                FArrayBox tmpfab(tmpbox, tmpcomp);
                Elixir eli = tmpfab.elixir();

                Array4<Real> fx = tmpfab.array(0);
                Array4<Real> fy = tmpfab.array(nmaxcomp);
                Array4<Real> fz = tmpfab.array(nmaxcomp*2);

                // velocity
                mol::compute_convective_fluxes(lev, bx, AMREX_SPACEDIM, fx, fy, fz, vel(lev).const_array(mfi),
                                          u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                                          vel.bcrec().data(),
                                          vel.bcrec_device().data(),geom);
                mol::compute_convective_rate(bx, AMREX_SPACEDIM, conv_u(lev).array(mfi), fx, fy, fz, geom[lev].InvCellSizeArray());

                // density
                if (!constant_density) {
                    mol::compute_convective_fluxes(lev, bx, 1, fx, fy, fz, den(lev).const_array(mfi),
                                              u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                                              den.bcrec().data(),
                                              den.bcrec_device().data(),geom);
                    mol::compute_convective_rate(bx, 1, conv_r(lev).array(mfi), fx, fy, fz, geom[lev].InvCellSizeArray());
                }

                // tracer
                if (advect_tracer) {
                    mol::compute_convective_fluxes(lev, bx, ntrac, fx, fy, fz, rhotrac,
                                              u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                                              tra.bcrec().data(),
                                              tra.bcrec_device().data(),geom);
                    mol::compute_convective_rate(bx, ntrac, conv_t(lev).array(mfi), fx, fy, fz, geom[lev].InvCellSizeArray());
                }

            }

        }
    }
}


