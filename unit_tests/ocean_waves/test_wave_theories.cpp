#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/relaxation_zones/stokes_waves_K.H"

namespace amr_wind_tests {

class WaveTheoriesTest : public MeshTest
{};

TEST_F(WaveTheoriesTest, StokesWaves)
{
    amrex::Real CoeffTol = 1e-4;

    amrex::Real wavenumber = 2.0;
    amrex::Real waterdepth = 0.376991;
    int StokesOrder = 5;
    amrex::Real c0, a11, a22, b22, c2, a31, a33, b31, a42;
    amrex::Real a44, b42, b44, c4, a51, a53, a55, b53, b55;

    amr_wind::ocean_waves::relaxation_zones::stokes_coefficients(
        StokesOrder, wavenumber, waterdepth, c0, a11, a22, b22, c2, a31, a33,
        b31, a42, a44, b42, b44, c4, a51, a53, a55, b53, b55);

    // Gold coefficient values taken from table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234

    const amrex::Real gold_A11 = 1.208490;
    const amrex::Real gold_A22 = 0.799840;
    const amrex::Real gold_A31 = -9.105340;
    const amrex::Real gold_A33 = 0.368275;
    const amrex::Real gold_A42 = -12.196150;
    const amrex::Real gold_A44 = 0.058723;
    const amrex::Real gold_A51 = 108.46831725;
    const amrex::Real gold_A53 = -6.941756;
    const amrex::Real gold_A55 = -0.074979;
    const amrex::Real gold_B22 = 2.502414;
    const amrex::Real gold_B31 = -5.731666;
    const amrex::Real gold_B42 = -32.407508;
    const amrex::Real gold_B44 = 14.033758;
    const amrex::Real gold_B53 = -103.44536875;
    const amrex::Real gold_B55 = 37.200027;
    const amrex::Real gold_C0 = 0.798448;
    const amrex::Real gold_C2 = 1.940215;
    const amrex::Real gold_C4 = -12.970403;

    EXPECT_NEAR(gold_A11, a11, CoeffTol);
    EXPECT_NEAR(gold_A22, a22, CoeffTol);
    EXPECT_NEAR(gold_A31, a31, CoeffTol);
    EXPECT_NEAR(gold_A33, a33, CoeffTol);
    EXPECT_NEAR(gold_A42, a42, CoeffTol);
    EXPECT_NEAR(gold_A44, a44, CoeffTol);
    EXPECT_NEAR(gold_A51, a51, CoeffTol);
    EXPECT_NEAR(gold_A53, a53, CoeffTol);
    EXPECT_NEAR(gold_A55, a55, CoeffTol);
    EXPECT_NEAR(gold_B22, b22, CoeffTol);
    EXPECT_NEAR(gold_B31, b31, CoeffTol);
    EXPECT_NEAR(gold_B42, b42, CoeffTol);
    EXPECT_NEAR(gold_B44, b44, CoeffTol);
    EXPECT_NEAR(gold_B53, b53, CoeffTol);
    EXPECT_NEAR(gold_B55, b55, CoeffTol);
    EXPECT_NEAR(gold_C0, c0, CoeffTol);
    EXPECT_NEAR(gold_C2, c2, CoeffTol);
    EXPECT_NEAR(gold_C4, c4, CoeffTol);
}

TEST_F(WaveTheoriesTest, StokesWavesFreeSurfaceProfile)
{
    amrex::Real Tol = 1e-4;

    amrex::Real g = 9.81;
    int StokesOrder = 5;
    // wavenumber k and waterdepth d chosen so that kd = 0.758932
    // to match value of column 3 from table 2 in
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    amrex::Real wavenumber = 2.0;
    amrex::Real waterdepth = 0.376991;
    amrex::Real wavelength = 2. * M_PI / wavenumber;
    amrex::Real waveheight = 0.1;
    amrex::Real zsl = 0.;
    amrex::Real x = 0.;
    amrex::Real z = -0.25;
    amrex::Real phase_offset = 0.;
    amrex::Real time = 0.;
    amrex::Real eta = 0.;
    amrex::Real u_w = 0.;
    amrex::Real v_w = 0.;
    amrex::Real w_w = 0.;

    amr_wind::ocean_waves::relaxation_zones::stokes_waves(
        StokesOrder, wavelength, waterdepth, waveheight, zsl, g, x, z, time,
        phase_offset, eta, u_w, v_w, w_w);

    // Coefficients values taken from column 3 of table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    amrex::Real B22 = 2.502414;
    amrex::Real B31 = -5.731666;
    amrex::Real B42 = -32.407508;
    amrex::Real B44 = 14.033758;
    amrex::Real B53 = -103.44536875;
    amrex::Real B55 = 37.200027;

    amrex::Real eps = wavenumber * waveheight / 2.;
    amrex::Real S = 2. * std::exp(2. * wavenumber * waterdepth) /
                    (std::exp(4. * wavenumber * waterdepth) + 1.);
    amrex::Real C = 1.0 - S;
    amrex::Real C0 = std::sqrt(std::tanh(wavenumber * waterdepth));
    amrex::Real C2 = C0 * (2 + 7 * std::pow(S, 2)) / (4 * std::pow(C, 2));
    amrex::Real C4 = C0 *
                     (4 + 32 * S - 116 * std::pow(S, 2) - 400 * std::pow(S, 3) -
                      71 * std::pow(S, 4) + 146 * std::pow(S, 5)) /
                     (32 * std::pow(C, 5));
    amrex::Real wave_speed =
        (C0 + std::pow(eps, 2) * C2 + std::pow(eps, 4) * C4) *
        std::sqrt(g / wavenumber);

    amrex::Real omega = wave_speed * wavenumber;
    amrex::Real phase = wavenumber * x - omega * time - phase_offset;

    // Check against Eq. (14) from Fenton 1985
    amrex::Real eta_theory =
        (eps * std::cos(phase) + std::pow(eps, 2) * B22 * std::cos(2. * phase) +
         std::pow(eps, 3) * B31 * (std::cos(phase) - std::cos(3. * phase)) +
         std::pow(eps, 4) *
             (B42 * std::cos(2. * phase) + B44 * std::cos(4. * phase)) +
         std::pow(eps, 5) *
             (-(B53 + B55) * std::cos(phase) + B53 * std::cos(3 * phase) +
              B55 * std::cos(5 * phase))) /
            wavenumber +
        zsl;

    EXPECT_NEAR(eta, eta_theory, Tol);

    // Re-evaluate with new set of coefficients
    // Deep-water limit (k*d->\infty)
    waterdepth = 100;
    waveheight = 0.16;
    wavenumber = 0.156;
    wavelength = 2. * M_PI / wavenumber;
    zsl = 0.;
    x = 4.;
    z = 0.;
    phase_offset = M_PI;
    time = 2.7;
    eta = 0.;
    u_w = 0.;
    v_w = 0.;
    w_w = 0.;

    amr_wind::ocean_waves::relaxation_zones::stokes_waves(
        StokesOrder, wavelength, waterdepth, waveheight, zsl, g, x, z, time,
        phase_offset, eta, u_w, v_w, w_w);

    // Coefficients values taken from column 1 of table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    B22 = 0.5;
    B31 = -0.375;
    B42 = 0.3333333;
    B44 = 0.3333333;
    B53 = 0.7734375;
    B55 = 0.3255208;

    eps = wavenumber * waveheight / 2.;
    // Coefficients computed analytically by taking limit kd -> \infty and
    // simplified accordingly
    // Note that in this limit S = 0 and C = 1 and thus they are omitted here
    C0 = 1.0;
    C2 = 0.5;
    C4 = 0.125;
    wave_speed = (C0 + std::pow(eps, 2) * C2 + std::pow(eps, 4) * C4) *
                 std::sqrt(g / wavenumber);

    omega = wave_speed * wavenumber;
    phase = wavenumber * x - omega * time - phase_offset;

    // Matches Eq. (18) from Fenton 1985
    eta_theory =
        (eps * std::cos(phase) + std::pow(eps, 2) * B22 * std::cos(2. * phase) +
         std::pow(eps, 3) * B31 * (std::cos(phase) - std::cos(3. * phase)) +
         std::pow(eps, 4) *
             (B42 * std::cos(2. * phase) + B44 * std::cos(4. * phase)) +
         std::pow(eps, 5) *
             (-(B53 + B55) * std::cos(phase) + B53 * std::cos(3 * phase) +
              B55 * std::cos(5 * phase))) /
            wavenumber +
        zsl;

    EXPECT_NEAR(eta, eta_theory, Tol);
}

TEST_F(WaveTheoriesTest, StokesWaveLength)
{
    // Values of wave_height, wave_period and water_depth taken from
    // https://www.sciencedirect.com/science/article/pii/S0029801817306066
    amrex::Real wave_height = 0.16;
    amrex::Real wave_period = 1.6;
    amrex::Real water_depth = 18.;
    int wave_order = 2;
    constexpr amrex::Real g = 9.81;
    constexpr amrex::Real tol_lambda = 1e-10;
    // To return computed wavelength from first guess of wavenumber k
    int iter_max = -1;

    // Check initial guess of wavenumber k for Newton iterations
    amrex::Real lambda =
        amr_wind::ocean_waves::relaxation_zones::stokes_wave_length(
            wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
            iter_max);

    const amrex::Real k_Newton = 2.0 * M_PI / lambda;

    // Compare with expected wavenumber from theory k = omega^2/g, where omega =
    // 2Pi/wave_period
    const amrex::Real k_theory =
        (2.0 * M_PI / wave_period) * (2.0 * M_PI / wave_period) / g;

    EXPECT_NEAR(k_Newton, k_theory, 1e-8);

    // Check wave theory
    wave_height = 0.2;
    wave_period = 2.2;
    water_depth = 5.0;
    wave_order = 5;
    iter_max = 20;

    lambda = amr_wind::ocean_waves::relaxation_zones::stokes_wave_length(
        wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
        iter_max);

    // Relation to check is Eq.(24) from course notes:
    // https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf
    amrex::Real k = 2.0 * M_PI / lambda;
    const amrex::Real RHS1 = 2.0 * M_PI / (wave_period * std::sqrt(g * k));

    amrex::Real S = 1.0 / std::cosh(2.0 * k * water_depth);
    amrex::Real C = 1.0 - S;
    amrex::Real eps = k * wave_height / 2.0;

    amrex::Real C0 = std::sqrt(std::tanh(k * water_depth));
    amrex::Real C2 = C0 * (2.0 + 7.0 * S * S) / (4.0 * C * C);
    const amrex::Real C4 =
        C0 *
        (4.0 + 32.0 * S - 116.0 * std::pow(S, 2) - 400.0 * std::pow(S, 3) -
         71.0 * std::pow(S, 4) + 146.0 * std::pow(S, 5)) /
        (32.0 * std::pow(C, 5));
    const amrex::Real LHS1 = C0 + std::pow(eps, 2) * C2 + std::pow(eps, 4) * C4;
    EXPECT_NEAR(LHS1, RHS1, 1e-8);

    // Reevaluate with a new set of conditions
    wave_height = 0.05;
    wave_period = 1.5;
    water_depth = 0.9;
    wave_order = 3;
    iter_max = 31;
    lambda = amr_wind::ocean_waves::relaxation_zones::stokes_wave_length(
        wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
        iter_max);

    k = 2.0 * M_PI / lambda;
    const amrex::Real RHS2 = 2.0 * M_PI / (wave_period * std::sqrt(g * k));

    S = 1.0 / std::cosh(2.0 * k * water_depth);
    C = 1.0 - S;
    eps = k * wave_height / 2.0;

    C0 = std::sqrt(std::tanh(k * water_depth));
    C2 = C0 * (2.0 + 7.0 * S * S) / (4.0 * C * C);
    const amrex::Real LHS2 = C0 + std::pow(eps, 2) * C2;
    EXPECT_NEAR(LHS2, RHS2, 1e-8);
}

} // namespace amr_wind_tests