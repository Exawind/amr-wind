#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/boundary_conditions/wall_models/MOSD.H"

namespace amr_wind_tests {

class MOSDTest : public MeshTest
{
protected:
    void populate_parameters() override { MeshTest::populate_parameters(); }
	
    amrex::Real get_u_loglaw(const amrex::Real utau) const
    {
        return utau * (1/0.4 * std::log(m_zref * utau / m_nu)  + 5);
    }

    // vertical location (z = dx)
    const amrex::Real m_zref = 6.2831853 / 128.0;
    // Molecular viscosity
    const amrex::Real m_nu = 1e-5;
    // Set test tolerance to convergence test in MOSD
    const amrex::Real m_tol = 1.0e-5;
};

TEST_F(MOSDTest, test_mosd)
{
    populate_parameters();
    amr_wind::MOSD md;

    md.amplitude = 0.05; 
    md.wavenumber = 4;
    md.omega = 0.8;
    md.time = 0;

    amrex::Real utau = 0.1;
    amrex::Real tau_wave_expected = 0.0;
    amrex::Real u_dx = get_u_loglaw(utau);
    amrex::Real v_dx = 0.0;
    amrex::Real x_c = 3.141592;  //For this location of the wave, the stress has to be zero, because sin(#pi)=zero//
    amrex::Real unit_nor = 0.0;
    auto tau_wave = md.get_dyn_tau(u_dx,v_dx,x_c,unit_nor);
    EXPECT_NEAR(tau_wave, tau_wave_expected, m_tol);

    
    tau_wave_expected = 0.0;
    x_c = 0.392699; //This give sin(pi/two) = one//
    //However, because this location is at the negative slope of the wave, the heaviside function will make the tau=zero
    EXPECT_NEAR(md.get_dyn_tau(u_dx,v_dx,x_c,unit_nor), tau_wave_expected, m_tol);

    tau_wave_expected = 0.02493;
    x_c = 1.0; //This give sin(four)= negative value
    //At this location we are in the possitive slope, therfore heaviside function will allow for stress
    EXPECT_NEAR(md.get_dyn_tau(u_dx,v_dx,x_c,unit_nor), tau_wave_expected, m_tol);

}
} // namespace amr_wind_tests
