#include "aw_test_utils/AmrexTest.H"

#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/utilities/trig_ops.H"

#include <string>

namespace amr_wind_tests {
namespace {

std::stringstream generate_txt_airfoil()
{
    using namespace ::amr_wind::utils;
    std::stringstream ss;
    ss << 6 << std::endl;
    ss << -180.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
    ss << -160.0 << " " << two_pi() * radians(20.0) << " " << 0.0 << " " << 0.0
       << std::endl;
    ss << -20.0 << " " << -two_pi() * radians(20.0) << " " << 0.0 << " " << 0.0
       << std::endl;
    ss << 20.0 << " " << two_pi() * radians(20.0) << " " << 0.0 << " " << 0.0
       << std::endl;
    ss << 160.0 << " " << -two_pi() * radians(20.0) << " " << 0.0 << " " << 0.0
       << std::endl;
    ss << 180.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;

    return ss;
}

std::stringstream generate_openfast_airfoil()
{
    std::stringstream ss;
    ss << "!........................................ " << std::endl;
    ss << "! Table of aerodynamics coefficients " << std::endl;
    ss << "        6   NumAlf            ! Number of data lines in the "
          "following table "
       << std::endl;
    ss << "!    Alpha      Cl      Cd        Cm  " << std::endl;
    ss << "!    (deg)      (-)     (-)       (-) " << std::endl;
    ss << "   -180.00    0.000   0.0407   0.0000 " << std::endl;
    ss << "   -175.00    0.223   0.0507   0.0937 " << std::endl;
    ss << "   -170.00    0.405   0.1055   0.1702 " << std::endl;
    ss << "   -160.00    0.658   0.2982   0.2819 " << std::endl;
    ss << "   -155.00    0.733   0.4121   0.3213 " << std::endl;
    ss << "   -150.00    0.778   0.5308   0.3520 " << std::endl;
    ss << "   -145.00    0.795   0.6503   0.3754 " << std::endl;
    ss << "   -140.00    0.787   0.7672   0.3926 " << std::endl;
    ss << "   -135.00    0.757   0.8785   0.4046 " << std::endl;

    return ss;
}

} // namespace

TEST(Airfoil, read_txt_file)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_txt_airfoil();

    auto af = AirfoilLoader::load_text_file(ss);
    EXPECT_EQ(af->num_entries(), 6);
    EXPECT_NEAR(af->aoa().front(), -1.0 * ::amr_wind::utils::pi(), 1.0e-12);
    EXPECT_NEAR(af->aoa().back(), 1.0 * ::amr_wind::utils::pi(), 1.0e-12);
}

TEST(Airfoil, read_openfast_file)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_openfast_airfoil();

    auto af = AirfoilLoader::load_openfast_airfoil(ss);
    EXPECT_EQ(af->num_entries(), 6);
    EXPECT_NEAR(af->aoa().front(), -1.0 * ::amr_wind::utils::pi(), 1.0e-12);
    EXPECT_NEAR(af->aoa().back(), ::amr_wind::utils::radians(-150.0), 1.0e-12);
}

TEST(Airfoil, airfoil_lookup)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_txt_airfoil();

    auto af = AirfoilLoader::load_text_file(ss);
    EXPECT_EQ(af->num_entries(), 6);

    amrex::Real cl, cd;
    amrex::Vector<amrex::Real> aoa_test{-15.0, -5.0, 0.0, 5.0, 10.0, 15.0};
    for (const auto& aoa_deg : aoa_test) {
        const amrex::Real aoa_rad = ::amr_wind::utils::radians(aoa_deg);
        (*af)(aoa_rad, cl, cd);
        EXPECT_NEAR(cl, ::amr_wind::utils::two_pi() * aoa_rad, 1.0e-3);
        EXPECT_NEAR(cd, 0.0, 1.0e-3);
    }
}

} // namespace amr_wind_tests
