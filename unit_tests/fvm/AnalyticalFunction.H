#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/MeshTest.H"

namespace amr_wind_tests {

class AnalyticalFunction
{
public:
    AnalyticalFunction(int degree, int ncomponents);

    ~AnalyticalFunction();

    amrex::Real
    phi_eval(amrex::Real x, amrex::Real y, amrex::Real z, int component) const;

    amrex::Real
    dphi_eval(amrex::Real x, amrex::Real y, amrex::Real z, int component) const;

    amrex::Real divergence(amrex::Real x, amrex::Real y, amrex::Real z) const;

    amrex::Real strainrate(amrex::Real x, amrex::Real y, amrex::Real z) const;

    amnrex::Real laplacian(amrex::Real x, amrex::Real y, amrex::Real z) const;

    const int m_degree;
    const int m_ncomponents;
    amrex::Vector<amrex::Real> m_coeffs;
};

} // namespace amr_wind_tests
