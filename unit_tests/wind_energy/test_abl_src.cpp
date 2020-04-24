#include "abl_test_utils.H"
#include "trig_ops.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "ABLForcing.H"
#include "CoriolisForcing.H"
#include "BoussinesqBuoyancy.H"

namespace amr_wind_tests {

TEST_F(ABLTest, abl_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;

    pp_utils::default_time_inputs();
    utils::populate_abl_params();

    amrex::Box bx{{0, 0, 0}, {2, 2, 2}};
    amrex::FArrayBox vel_src(bx, AMREX_SPACEDIM);

    amr_wind::SimTime time;
    time.parse_parameters();

    amr_wind::ABLForcingOld abl_forcing(time);

    // During initialization ensure that the source terms are zero
    vel_src.setVal<amrex::RunOn::Device>(0.0);
    abl_forcing(bx, vel_src.array());
    for (int i=0; i < AMREX_SPACEDIM; ++i) {
        const auto min_src = vel_src.min<amrex::RunOn::Device>(i);
        const auto max_src = vel_src.max<amrex::RunOn::Device>(i);
        // Ensure that the source term is zero for initialization
        EXPECT_NEAR(min_src, 0.0, tol);
        // Ensure that the source term is constant throughout the domain
        EXPECT_NEAR(min_src, max_src, tol);
    }

    // Mimic source term at later timesteps
    {
        time.new_timestep();
        time.set_current_cfl(2.0);
        EXPECT_NEAR(time.deltaT(), 0.1, tol);

        vel_src.setVal<amrex::RunOn::Device>(0.0);
        abl_forcing.set_mean_velocities(10.0, 5.0);
        abl_forcing(bx, vel_src.array());

        // Targets are U = (20.0, 10.0, 0.0) set in initial conditions
        // Means (set above) V = (10.0, 5.0, 0.0)
        // deltaT (set above) dt = 0.1
        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{{100.0, 50.0, 0.0}};
        for (int i=0; i < AMREX_SPACEDIM; ++i) {
            const auto min_src = vel_src.min<amrex::RunOn::Device>(i);
            const auto max_src = vel_src.max<amrex::RunOn::Device>(i);
            EXPECT_NEAR(min_src, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_src, max_src, tol);
        }
    }
}

TEST_F(ABLTest, coriolis_const_vel)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));
    // Initialize a random value for the velocity component
    const amrex::Real vel_comp = 10.0 + 5.0 * (amrex::Random() - 0.5);

    // Initialize parameters
    utils::populate_abl_params();

    // Create velocity and source terms fields
    amrex::Box bx{{0, 0, 0}, {2, 2, 2}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox vel_src(bx, AMREX_SPACEDIM);

    amr_wind::CoriolisForcingOld coriolis;


    // Velocity in x-direction test
    {
        amrex::Real golds[AMREX_SPACEDIM] = {0.0, -corfac * latfac * vel_comp,
                                             corfac * latfac * vel_comp};
        // Reset fields
        velocity.setVal<amrex::RunOn::Device>(0.0);
        vel_src.setVal<amrex::RunOn::Device>(0.0);
        // set x component
        velocity.setVal<amrex::RunOn::Device>(vel_comp, 0);

        coriolis(bx, velocity.array(), vel_src.array());

        for (int i=0; i < AMREX_SPACEDIM; ++i) {
            const auto min_src = vel_src.min<amrex::RunOn::Device>(i);
            const auto max_src = vel_src.max<amrex::RunOn::Device>(i);
            EXPECT_NEAR(min_src, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_src, max_src, tol);
        }
    }

    // Velocity in y direction test
    {
        amrex::Real golds[AMREX_SPACEDIM] = {corfac * latfac * vel_comp, 0.0, 0.0};
        // Reset fields
        velocity.setVal<amrex::RunOn::Device>(0.0);
        vel_src.setVal<amrex::RunOn::Device>(0.0);
        // Set y component
        velocity.setVal<amrex::RunOn::Device>(vel_comp, 1);

        coriolis(bx, velocity.array(), vel_src.array());

        for (int i=0; i < AMREX_SPACEDIM; ++i) {
            const auto min_src = vel_src.min<amrex::RunOn::Device>(i);
            const auto max_src = vel_src.max<amrex::RunOn::Device>(i);
            EXPECT_NEAR(min_src, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_src, max_src, tol);
        }
    }
}

namespace {

void cor_height_init_vel_field(
    amrex::Box& bx,
    amrex::FArrayBox& velocity)
{
    // Set y velocity as a function of height with (dx = 1.0)
    auto vel = velocity.array();
    amrex::ParallelFor(bx, [vel] AMREX_GPU_DEVICE(int i, int j, int k) {
        vel(i, j, k, 1) = static_cast<amrex::Real>(k);
    });
}

}

TEST_F(ABLTest, coriolis_height_variation)
{
    constexpr int kdim = 16;
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));

    // Initialize parameters
    utils::populate_abl_params();

    // Create velocity and source terms fields
    amrex::Box bx{{0, 0, 0}, {2, 2, kdim}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox vel_src(bx, AMREX_SPACEDIM);

    amr_wind::CoriolisForcingOld coriolis;

    velocity.setVal<amrex::RunOn::Device>(0.0);
    vel_src.setVal<amrex::RunOn::Device>(0.0);
    cor_height_init_vel_field(bx, velocity);

    coriolis(bx, velocity.array(), vel_src.array());

    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(0), corfac * latfac * kdim, tol);

    for (int i=1; i < AMREX_SPACEDIM; ++i) {
        const auto min_src = vel_src.min<amrex::RunOn::Device>(i);
        const auto max_src = vel_src.max<amrex::RunOn::Device>(i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(min_src, max_src, tol);
    }
}

namespace {

void init_abl_temperature_field(int kdim,
                           amrex::Box& bx,
                           amrex::FArrayBox& tracer)
{
    // Set tracer as a function of height with (dx = 1.0)
    auto trac = tracer.array();
    const amrex::Real dz = 1000.0/((amrex::Real) kdim+1);
    
    amrex::ParallelFor(bx, [dz,trac] AMREX_GPU_DEVICE(int i, int j, int k) {

        const amrex::Real z = (k+0.5)*dz;
        
        // potential temperature profile
        if(z < 650.0){
            trac(i, j, k, 0) = 300.0;
        } else if(z < 750.0){
            trac(i, j, k, 0) = 300.0 + (z-650.0)/(750.0-650.0)*8.0;
        } else {
            trac(i, j, k, 0) = 308.0;
        }
        
    });
}

}

TEST_F(ABLTest, boussinesq)
{
    constexpr int kdim = 16;
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_abl_params();

    // Create velocity and source terms fields
    amrex::Box bx{{0, 0, 0}, {2, 2, kdim}};
    amrex::FArrayBox temperature(bx, 1);
    amrex::FArrayBox vel_src(bx, AMREX_SPACEDIM);

    amr_wind::BoussinesqBuoyancyOld bb;
    
    temperature.setVal<amrex::RunOn::Device>(0.0);
    vel_src.setVal<amrex::RunOn::Device>(0.0);
    init_abl_temperature_field(kdim, bx, temperature);

    bb(bx, temperature.array(), vel_src.array());

    // should be no forcing in x and y directions
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(1), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(1), 0.0, tol);

//    f = beta * (T0 - T)*g
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(2), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(2), -9.81*(300.0-308.0)/300.0, tol);

}

#if 0

namespace {

void init_density_field(int kdim,
                           amrex::Box& bx,
                           amrex::FArrayBox& density)
{
    // Set density as a function of height with (dz = 1.0)
    auto den = density.array();
    const amrex::Real dz = 1.0/((amrex::Real) kdim+1);

    amrex::ParallelFor(bx, [dz,den] AMREX_GPU_DEVICE(int i, int j, int k) {

        const amrex::Real z = (k+0.5)*dz;

        if(z < 0.3){
            den(i, j, k) = 0.5;
        } else if(z > 0.7){
            den(i, j, k) = 2.0;
        } else {
            den(i, j, k) = 1.0;
        }

    });
}

}

TEST_F(ABLTest, densitybuoyancy)
{
    constexpr int kdim = 16;
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_abl_params();

    // Create velocity and source terms fields
    amrex::Box bx{{0, 0, 0}, {2, 2, kdim}};
    amrex::FArrayBox density(bx, 1);
    amrex::FArrayBox vel_src(bx, AMREX_SPACEDIM);
    amr_wind::DensityBuoyancy db;

    density.setVal<amrex::RunOn::Device>(1.0);
    vel_src.setVal<amrex::RunOn::Device>(0.0);

    db(bx, density.array(), vel_src.array());

    // first test constant density field

    // should be no forcing in x and y directions
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(1), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(1), 0.0, tol);

    // f = (1-rho_0/rho)*g
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(2), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(2), 0.0, tol);

    init_density_field(kdim, bx, density);
    vel_src.setVal<amrex::RunOn::Device>(0.0);
    db(bx, density.array(), vel_src.array());

    // should be no forcing in x and y directions
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(0), 0.0, tol);
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(1), 0.0, tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(1), 0.0, tol);

    // f = g*(1-rho_0/rho)
    EXPECT_NEAR(vel_src.min<amrex::RunOn::Device>(2), -9.81*(1.0-1.0/2.0), tol);
    EXPECT_NEAR(vel_src.max<amrex::RunOn::Device>(2), -9.81*(1.0-1.0/0.5), tol);
}
#endif

} // namespace amr_wind_tests
