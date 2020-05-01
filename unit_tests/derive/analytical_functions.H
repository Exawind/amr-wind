/** \file analytical_functions.cpp
 *
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"

namespace amr_wind_tests {

//! Create unique namespace for this test fixture
class AnalyticalFunctions  
{

public:

    AnalyticalFunctions(int n, amrex::Box);

    ~AnalyticalFunctions();

    amrex::FArrayBox scalar;
    amrex::FArrayBox scalargrad;

    amrex::Real dx,dy,dz;

};

AnalyticalFunctions::AnalyticalFunctions(int n,amrex::Box bx)
: scalar(bx,1), scalargrad(bx,AMREX_SPACEDIM)
{
    // Define scalar fields and gradient
    dx = 0.1; dy = 0.2 + 0.01*amrex::Random(); dz = 0.3;
    
    auto phi = scalar.array();
    auto gradphi = scalargrad.array();
    
    // Interior stencil should be exact for quadratic polynomials
    for(int i = 0; i <= n; ++i){
        for(int j = 0; j <= n; ++j){
            for(int k = 0; k <= n; ++k){
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5 + j*dy;
                amrex::Real z = 5.2 + k*dz;
                
                phi(i,j,k) = 5.0*x*y*z + 4.0*x*x + 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z + 9.4*x*y + 3.4*x+ 2.0*x + 3.14;
                gradphi(i,j,k,0) = 5.0*y*z + 8.0*x + 3.8*z + 9.4*y + 3.4 + 2.0; 
                gradphi(i,j,k,1) = 5.0*x*z + 6.0*y + 1.3*z + 9.4*x; 
                gradphi(i,j,k,2) = 5.0*x*y -4.6*z + 1.3*y + 3.8*x;
            }
        }
    }
}

AnalyticalFunctions::~AnalyticalFunctions()=default;

}
