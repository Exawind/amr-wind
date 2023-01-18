#include "abl_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "amr-wind/equation_systems/icns/source_terms/DensityBuoyancy.H"

namespace amr_wind_tests {

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(ABLMeshTest, abl_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    utils::populate_abl_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;

    amr_wind::pde::icns::ABLForcing abl_forcing(sim());

    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, 0.0, tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Mimic source term at later timesteps
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(2.0, 0.0, 0.0);
        EXPECT_NEAR(time.deltaT(), 0.1, tol);

        src_term.setVal(0.0);
        abl_forcing.set_mean_velocities(10.0, 5.0);
        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        // Targets are U = (20.0, 10.0, 0.0) set in initial conditions
        // Means (set above) V = (10.0, 5.0, 0.0)
        // deltaT (set above) dt = 0.1
        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {100.0, 50.0, 0.0}};
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }
}

TEST_F(ABLMeshTest, body_force)
{
    constexpr amrex::Real tol = 1.0e-12;
    utils::populate_abl_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;

    amr_wind::pde::icns::BodyForce body_force(sim());

    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        body_force(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    const auto valx = utils::field_max(src_term, 0);
    const auto valy = utils::field_max(src_term, 1);
    const auto valz = utils::field_max(src_term, 2);
    EXPECT_NEAR(valx, 1.0, tol);
    EXPECT_NEAR(valy, 2.0, tol);
    EXPECT_NEAR(valz, 3.0, tol);

    // Mimic source term at later timesteps
    {
        auto& time = sim().time();
        time.current_time() = 0.1;
        src_term.setVal(0.0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            body_force(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {static_cast<amrex::Real>(1.0 * std::cos(0.1)),
             static_cast<amrex::Real>(2.0 * std::cos(0.1)),
             static_cast<amrex::Real>(3.0 * std::cos(0.1))}};
        const auto valx2 = utils::field_max(src_term, 0);
        const auto valy2 = utils::field_max(src_term, 1);
        const auto valz2 = utils::field_max(src_term, 2);
        EXPECT_NEAR(valx2, golds[0], tol);
        EXPECT_NEAR(valy2, golds[1], tol);
        EXPECT_NEAR(valz2, golds[2], tol);
    }
}

TEST_F(ABLMeshTest, geostrophic_two_component_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));

    utils::populate_abl_params();

    amrex::ParmParse pp("CoriolisForcing");
    pp.add("latitude", 45.0);

    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& density = sim().repo().get_field("density");

    amr_wind::pde::icns::GeostrophicForcing geostrophic_forcing(sim());

    // Two component forcing
    {
        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {-corfac * 6.0 * latfac, corfac * 10.0 * latfac, 0.0}};

        density.setVal(1.0);
        src_term.setVal(0.0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            geostrophic_forcing(
                lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }
}

TEST_F(ABLMeshTest, geostrophic_three_component_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));

    utils::populate_abl_params();

    {
        amrex::ParmParse pp("ABL");
        pp.add("three_ComponentForcing", true);
    }

    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 45.0);
    }

    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& density = sim().repo().get_field("density");

    amr_wind::pde::icns::GeostrophicForcing geostrophic_forcing(sim());

    // Three component forcing
    {
        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {-corfac * 6.0 * latfac, +corfac * 10.0 * latfac,
             -corfac * 10.0 * latfac}};

        density.setVal(1.0);
        src_term.setVal(0.0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            geostrophic_forcing(
                lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }
}

TEST_F(ABLMeshTest, coriolis_two_component_const_vel)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));
    // Initialize a random value for the velocity component
    const amrex::Real vel_comp = 10.0 + 5.0 * (amrex::Random() - 0.5);
    amrex::ParmParse pp("CoriolisForcing");
    pp.add("latitude", 45.0);

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto fields = ICNSFields(sim())(sim().time());
    auto& vel = fields.field;
    auto& src_term = fields.src_term;
    amr_wind::pde::icns::CoriolisForcing coriolis(sim());

    // Velocity in x-direction test two component forcing
    {
        amrex::Real golds[AMREX_SPACEDIM] = {
            0.0, -corfac * latfac * vel_comp, 0.0};
        vel.setVal(0.0);
        src_term.setVal(0.0);
        vel.setVal(vel_comp, 0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            coriolis(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }

    // Velocity in y-direction test two component forcing
    {
        amrex::Real golds[AMREX_SPACEDIM] = {
            corfac * latfac * vel_comp, 0.0, 0.0};
        vel.setVal(0.0);
        src_term.setVal(0.0);
        vel.setVal(vel_comp, 1);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            coriolis(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }
}

TEST_F(ABLMeshTest, coriolis_three_component_const_vel)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));
    // Initialize a random value for the velocity component
    const amrex::Real vel_comp = 10.0 + 5.0 * (amrex::Random() - 0.5);

    {
        amrex::ParmParse pp("ABL");
        pp.add("three_ComponentForcing", true);
    }

    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 45.0);
    }

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto fields = ICNSFields(sim())(sim().time());
    auto& vel = fields.field;
    auto& src_term = fields.src_term;
    amr_wind::pde::icns::CoriolisForcing coriolis(sim());

    // Velocity in x-direction test three component forcing
    {
        amrex::Real golds[AMREX_SPACEDIM] = {
            0.0, -corfac * latfac * vel_comp, corfac * latfac * vel_comp};
        vel.setVal(0.0);
        src_term.setVal(0.0);
        vel.setVal(vel_comp, 0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            coriolis(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, golds[i], tol);
            // Ensure that the source term is constant throughout the domain
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }
}

namespace {

void cor_height_init_vel_field(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& vel)
{
    // Set y velocity as a function of height with (dx = 1.0)
    amrex::ParallelFor(bx, [vel] AMREX_GPU_DEVICE(int i, int j, int k) {
        vel(i, j, k, 1) = static_cast<amrex::Real>(k);
    });
}

} // namespace

TEST_F(ABLMeshTest, coriolis_height_variation)
{
    constexpr int kdim = 7;
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto fields = ICNSFields(sim())(sim().time());
    auto& velocity = fields.field;
    auto& vel_src = fields.src_term;
    amr_wind::pde::icns::CoriolisForcing coriolis(sim());

    velocity.setVal(0.0);
    vel_src.setVal(0.0);

    run_algorithm(velocity, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& vel_arr = velocity(lev).array(mfi);
        const auto& vel_src_arr = vel_src(lev).array(mfi);
        cor_height_init_vel_field(bx, vel_arr);
        coriolis(lev, mfi, bx, amr_wind::FieldState::New, vel_src_arr);
    });

    EXPECT_NEAR(utils::field_min(vel_src, 0), 0.0, tol);
    EXPECT_NEAR(utils::field_max(vel_src, 0), corfac * latfac * kdim, tol);

    for (int i = 1; i < AMREX_SPACEDIM; ++i) {
        const auto min_src = utils::field_min(vel_src, i);
        const auto max_src = utils::field_max(vel_src, i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(min_src, max_src, tol);
    }
}

namespace {

void init_abl_temperature_field(
    int kdim, const amrex::Box& bx, const amrex::Array4<amrex::Real>& trac)
{
    // Set tracer as a function of height with (dx = 1.0)
    const amrex::Real dz = 1000.0 / ((amrex::Real)kdim + 1);

    amrex::ParallelFor(bx, [dz, trac] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real z = (k + 0.5) * dz;

        // potential temperature profile
        if (z < 650.0) {
            trac(i, j, k, 0) = 300.0;
        } else if (z < 750.0) {
            trac(i, j, k, 0) = 300.0 + (z - 650.0) / (750.0 - 650.0) * 8.0;
        } else {
            trac(i, j, k, 0) = 308.0;
        }
    });
}

} // namespace

TEST_F(ABLMeshTest, boussinesq)
{
    constexpr int kdim = 7;
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");
    sim().init_physics();

    amr_wind::pde::icns::BoussinesqBuoyancy bb(sim());

    auto& src_term = pde_mgr.icns().fields().src_term;

    auto& temperature =
        sim().repo().get_field("temperature", amr_wind::FieldState::Old);

    src_term.setVal(0.0);

    run_algorithm(temperature, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& temp_arr = temperature(lev).array(mfi);
        const auto& src_arr = src_term(lev).array(mfi);

        init_abl_temperature_field(kdim, bx, temp_arr);
        bb(lev, mfi, bx, amr_wind::FieldState::Old, src_arr);
    });

    // should be no forcing in x and y directions
    for (int i = 0; i < 2; ++i) {
        const auto min_src = utils::field_min(src_term, i);
        const auto max_src = utils::field_max(src_term, i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(max_src, 0.0, tol);
    }

    // f = beta * (T0 - T)*g
    EXPECT_NEAR(utils::field_min(src_term, 2), 0.0, tol);
    EXPECT_NEAR(
        utils::field_max(src_term, 2), -9.81 * (300.0 - 308.0) / 300.0, tol);
}

TEST_F(ABLMeshTest, boussinesq_nph)
{
    constexpr int kdim = 7;
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");
    sim().init_physics();

    amr_wind::pde::icns::BoussinesqBuoyancy bb(sim());

    auto& src_term = pde_mgr.icns().fields().src_term;

    auto& temperature =
        sim().repo().get_field("temperature", amr_wind::FieldState::NPH);

    src_term.setVal(0.0);

    run_algorithm(temperature, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& temp_arr = temperature(lev).array(mfi);
        const auto& src_arr = src_term(lev).array(mfi);

        init_abl_temperature_field(kdim, bx, temp_arr);
        bb(lev, mfi, bx, amr_wind::FieldState::NPH, src_arr);
    });

    // should be no forcing in x and y directions
    for (int i = 0; i < 2; ++i) {
        const auto min_src = utils::field_min(src_term, i);
        const auto max_src = utils::field_max(src_term, i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(max_src, 0.0, tol);
    }

    // f = beta * (T0 - T)*g
    EXPECT_NEAR(utils::field_min(src_term, 2), 0.0, tol);
    EXPECT_NEAR(
        utils::field_max(src_term, 2), -9.81 * (300.0 - 308.0) / 300.0, tol);
}

namespace {

void init_density_field(
    int kdim, const amrex::Box& bx, const amrex::Array4<amrex::Real>& den)
{
    // Set density as a function of height with (dz = 1.0)
    const amrex::Real dz = 1.0 / ((amrex::Real)kdim + 1);

    amrex::ParallelFor(bx, [dz, den] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real z = (k + 0.5) * dz;

        if (z < 0.3) {
            den(i, j, k) = 0.5;
        } else if (z > 0.7) {
            den(i, j, k) = 2.0;
        } else {
            den(i, j, k) = 1.0;
        }
    });
}

} // namespace

TEST_F(ABLMeshTest, densitybuoyancy)
{
    constexpr int kdim = 7;
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_abl_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Density");
    pde_mgr.register_transport_pde("Temperature");
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& density =
        sim().repo().get_field("density", amr_wind::FieldState::Old);

    amr_wind::pde::icns::DensityBuoyancy db(sim());

    density.setVal(1.0);
    src_term.setVal(0.0);

    run_algorithm(density, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& src_arr = src_term(lev).array(mfi);
        db(lev, mfi, bx, amr_wind::FieldState::Old, src_arr);
    });

    // should be no forcing for constant density
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_src = utils::field_min(src_term, i);
        const auto max_src = utils::field_max(src_term, i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(max_src, 0.0, tol);
    }

    run_algorithm(density, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& src_arr = src_term(lev).array(mfi);
        const auto& den_arr = density(lev).array(mfi);
        init_density_field(kdim, bx, den_arr);
        db(lev, mfi, bx, amr_wind::FieldState::Old, src_arr);
    });

    // should be no forcing in x and y directions
    for (int i = 0; i < 2; ++i) {
        const auto min_src = utils::field_min(src_term, i);
        const auto max_src = utils::field_max(src_term, i);
        EXPECT_NEAR(min_src, 0.0, tol);
        EXPECT_NEAR(max_src, 0.0, tol);
    }

    // f = g*(1-rho_0/rho)
    EXPECT_NEAR(utils::field_min(src_term, 2), -9.81 * (1.0 - 1.0 / 2.0), tol);
    EXPECT_NEAR(utils::field_max(src_term, 2), -9.81 * (1.0 - 1.0 / 0.5), tol);
}

} // namespace amr_wind_tests