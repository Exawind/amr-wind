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
#include "amr-wind/equation_systems/icns/source_terms/HurricaneForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/RayleighDamping.H"

#include "amr-wind/equation_systems/temperature/source_terms/HurricaneTempForcing.H"

namespace amr_wind_tests {

namespace {
amrex::Real get_val_at_height(
    amr_wind::Field& field,
    const int lev,
    const int comp,
    const amrex::Real ploz,
    const amrex::Real dz,
    const amrex::Real height)
{
    amrex::Real error_total = 0;

    error_total += amrex::ReduceSum(
        field(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                const amrex::Real z = ploz + (0.5 + k) * dz;
                // Check if current cell is closest to desired height
                if (std::abs(z - height) < 0.5 * dz) {
                    // Add field value to output
                    error += f_arr(i, j, k, comp);
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}
} // namespace

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(ABLMeshTest, abl_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();
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
        EXPECT_NEAR(time.delta_t(), 0.1, tol);

        src_term.setVal(0.0);
        abl_forcing.set_mean_velocities(10.0, 5.0);
        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        // Targets are U = (20.0, 10.0, 0.0) set in initial conditions
        // Means (set above) V = (10.0, 5.0, 0.0)
        // delta_t (set above) dt = 0.1
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
    populate_parameters();
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
        time.new_timestep();
        time.current_time() = 0.05;
        src_term.setVal(0.0);

        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            body_force(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {static_cast<amrex::Real>(
                 1.0 * std::cos(0.5 * (0.05 + time.new_time()))),
             static_cast<amrex::Real>(
                 2.0 * std::cos(0.5 * (0.05 + time.new_time()))),
             static_cast<amrex::Real>(
                 3.0 * std::cos(0.5 * (0.05 + time.new_time())))}};
        const auto valx2 = utils::field_max(src_term, 0);
        const auto valy2 = utils::field_max(src_term, 1);
        const auto valz2 = utils::field_max(src_term, 2);
        EXPECT_NEAR(valx2, golds[0], tol);
        EXPECT_NEAR(valy2, golds[1], tol);
        EXPECT_NEAR(valz2, golds[2], tol);
    }
}

TEST_F(ABLMeshTest, geostrophic_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();

    amrex::ParmParse pp("CoriolisForcing");
    pp.add("latitude", 54.0);

    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0);

    amr_wind::pde::icns::GeostrophicForcing geostrophic_forcing(sim());
    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        geostrophic_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    constexpr amrex::Real corfac =
        2.0 * amr_wind::utils::two_pi() / 86164.091 * 0.80901699437;
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
        {-corfac * 6.0, corfac * 10.0, 0.0}};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, golds[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }
}

TEST_F(ABLMeshTest, rayleigh_damping)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& velocity = sim().repo().get_field("velocity");
    velocity.setVal({{0.0, 0.0, 0.0}});
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0);

    amr_wind::pde::icns::RayleighDamping rayleigh_damping(sim());
    src_term.setVal(0.0);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        rayleigh_damping(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Domain is from 0 to 1000 in z, with 64 cells
    // Damping where coefficient is 1 is 50 long
    // Intermediate damping region (0 to 1) is 250 long
    const amrex::Real dz = sim().mesh().Geom(0).CellSizeArray()[2];
    const amrex::Real ploz = sim().mesh().Geom(0).ProbLoArray()[2];
    const amrex::Real phiz = sim().mesh().Geom(0).ProbHiArray()[2];
    // Testing locations:
    // 0: below top of domain
    // 1: above bottom of complete damping zone
    // 2: below top of graded damping zone
    // 3: above bottom of graded damping zone
    // 4: below bottom of graded damping (no damping)
    const amrex::Array<amrex::Real, 5> test_heights{
        phiz - 0.1, phiz - 50 + 0.5 * dz,
        phiz - dz * (0.5 + std::ceil(50 / dz)),
        phiz - dz * (-0.5 + std::ceil(250 / dz)), phiz - 250 - 0.5 * dz};
    // Expected values of coeff for each location
    const amrex::Array<amrex::Real, 5> golds{
        1., 1.,
        0.5 * std::cos(M_PI * (1000 - 50 - test_heights[2]) / 200) + 0.5,
        0.5 * std::cos(M_PI * (1000 - 50 - test_heights[3]) / 200) + 0.5, 0.0};

    // Get damping values from src term
    amrex::Array<amrex::Real, 5> src_x_vals;
    amrex::Array<amrex::Real, 5> src_y_vals;
    amrex::Array<amrex::Real, 5> src_z_vals;
    // Multiply to remove tau factor and divide by nx*ny cells because of
    // sum
    const amrex::Real tau = 40.;
    const int nx = 8;
    const int ny = 8;
    for (int n = 0; n < 5; ++n) {
        src_x_vals[n] =
            get_val_at_height(src_term, 0, 0, ploz, dz, test_heights[n]) * tau /
            nx / ny;
        src_y_vals[n] =
            get_val_at_height(src_term, 0, 1, ploz, dz, test_heights[n]) * tau /
            nx / ny;
        src_z_vals[n] =
            get_val_at_height(src_term, 0, 2, ploz, dz, test_heights[n]) * tau /
            nx / ny;
    }

    // Check each src value against expectations
    // Reference velocity is (12, 1, -3)
    // Damping is turned off in the y direction
    for (int n = 0; n < 5; ++n) {
        EXPECT_NEAR(src_x_vals[n], 12. * golds[n], tol);
        EXPECT_NEAR(src_y_vals[n], 1. * golds[n] * 0., tol);
        EXPECT_NEAR(src_z_vals[n], -3. * golds[n], tol);
    }
    // Check intermediate points against manual expectations (near 1, near
    // 0)
    EXPECT_NEAR(src_x_vals[2], 12. * 1.0, 12. * 0.05);
    EXPECT_NEAR(src_y_vals[2], 1. * 1.0 * 0.0, tol);
    EXPECT_NEAR(src_z_vals[2], -3. * 1.0, 3. * 0.05);
    EXPECT_NEAR(src_x_vals[3], 12. * 0.0, 12. * 0.005);
    EXPECT_NEAR(src_y_vals[3], 1. * 0.0 * 0.0, tol);
    EXPECT_NEAR(src_z_vals[3], -3. * 0.0, 3. * 0.005);
}

TEST_F(ABLMeshTest, hurricane_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();
    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 90.0);
    }
    {
        amrex::ParmParse pp("HurricaneTempForcing");
        pp.add("radial_decay", 0.001);
    }
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    auto& src_term_icns = pde_mgr.icns().fields().src_term;

    auto& temp_eqn = pde_mgr("Temperature-Godunov");
    auto& src_term_temp = temp_eqn.fields().src_term;

    auto& velocity = sim().repo().get_field("velocity");
    velocity.setVal({{0.0, 0.0, 0.0}});
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0);

    // Calculate Hurricane Forcing to the momentum equation
    amr_wind::pde::icns::HurricaneForcing hurricane_forcing(sim());
    src_term_icns.setVal(0.0);
    run_algorithm(src_term_icns, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term_icns(lev).array(mfi);

        hurricane_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Calculate Hurricane Forcing to the temperature equation
    amr_wind::pde::temperature::HurricaneTempForcing hurricane_temp_forcing(
        sim());
    src_term_temp.setVal(0.0);
    run_algorithm(src_term_temp, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term_temp(lev).array(mfi);

        hurricane_temp_forcing(
            lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Test the momentum hurricane forcing
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86400.0;
    const amrex::Real dz = sim().mesh().Geom(0).CellSizeArray()[2];
    const amrex::Real ratio_top = (18000. - (1000. - 0.5 * dz)) / 18000.;
    const amrex::Real ratio_bottom = (18000. - 0.5 * dz) / 18000.;
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds_max{
        {-corfac * 40.0 * ratio_top -
             40.0 * ratio_top * 40.0 * ratio_top / 40000.0,
         0.0, 0.0}};
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds_min{
        {-corfac * ratio_bottom * 40 -
             ratio_bottom * 40.0 * ratio_bottom * 40.0 / 40000.0,
         0.0, 0.0}};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto max_val = utils::field_max(src_term_icns, i);
        const auto min_val = utils::field_min(src_term_icns, i);
        EXPECT_NEAR(max_val, golds_max[i], tol);
        EXPECT_NEAR(min_val, golds_min[i], tol);
    }

    // Test the temperature equation hurricane forcing
    const amrex::Real gold = 0.0;
    const auto max_val = utils::field_max(src_term_temp);
    EXPECT_NEAR(max_val, gold, tol);
}

TEST_F(ABLMeshTest, coriolis_const_vel)
{
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86164.091;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));
    // Initialize a random value for the velocity component
    const amrex::Real vel_comp = 10.0 + 5.0 * (amrex::Random() - 0.5);

    // Initialize parameters
    populate_parameters();
    initialize_mesh();

    auto fields = ICNSFields(sim())(sim().time());
    auto& vel = fields.field;
    auto& src_term = fields.src_term;
    amr_wind::pde::icns::CoriolisForcing coriolis(sim());

    // Velocity in x-direction test
    {
        const amrex::RealArray golds = {
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

    // Velocity in y-direction test
    {
        const amrex::RealArray golds = {corfac * latfac * vel_comp, 0.0, 0.0};
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
    // ABL unit test mesh has 64 cells in z
    constexpr int kdim = 63;
    constexpr amrex::Real tol = 1.0e-12;
    constexpr amrex::Real corfac = 2.0 * amr_wind::utils::two_pi() / 86164.091;
    // Latitude is set to 45 degrees in the input file so sinphi = cosphi
    const amrex::Real latfac = std::sin(amr_wind::utils::radians(45.0));

    // Initialize parameters
    populate_parameters();
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
    populate_parameters();
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
    populate_parameters();
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
    populate_parameters();
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
