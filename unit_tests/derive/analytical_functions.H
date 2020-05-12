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

    void linear();
    
    void quadratic();

    int ncells_;
    amrex::Real dx_,dy_,dz_;
    amrex::FArrayBox scalar_;
    amrex::FArrayBox scalargrad_;

};

AnalyticalFunctions::AnalyticalFunctions(int n,amrex::Box bx)
  : ncells_(n),
    dx_(0.1),
    dy_(0.2 + 0.01*amrex::Random()),
    dz_(0.3),
    scalar_(bx,1),
    scalargrad_(bx,AMREX_SPACEDIM)
{}

void AnalyticalFunctions::quadratic()
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    
    // Interior stencil should be exact for quadratic polynomials
    for(int i = 0; i <= ncells_; ++i){
        for(int j = 0; j <= ncells_; ++j){
            for(int k = 0; k <= ncells_; ++k){
                amrex::Real x = 2.0 + i*dx_;
                amrex::Real y = 2.5 + j*dy_;
                amrex::Real z = 5.2 + k*dz_;
                
                phi(i,j,k) = 5.0*x*y*z + 4.0*x*x + 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z + 9.4*x*y + 3.4*x+ 2.0*x + 3.14;
                gradphi(i,j,k,0) = 5.0*y*z + 8.0*x + 3.8*z + 9.4*y + 3.4 + 2.0; 
                gradphi(i,j,k,1) = 5.0*x*z + 6.0*y + 1.3*z + 9.4*x; 
                gradphi(i,j,k,2) = 5.0*x*y -4.6*z + 1.3*y + 3.8*x;
            }
        }
    }
}

void AnalyticalFunctions::linear()
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    
    // Boundary stencil should be exact for linear polynomials
    for(int i = 0; i <= ncells_; ++i){
        for(int j = 0; j <= ncells_; ++j){
            for(int k = 0; k <= ncells_; ++k){
                amrex::Real x = 2.0 + i*dx_;
                amrex::Real y = 2.5 + j*dy_;
                amrex::Real z = 5.2 + k*dz_;
                
                phi(i,j,k) = 3.4*x+ 2.0*y + 7.8*z + 3.14;
                gradphi(i,j,k,0) = 3.4;
                gradphi(i,j,k,1) = 2.0;
                gradphi(i,j,k,2) = 7.8;
            }
        }
    }
    
}
AnalyticalFunctions::~AnalyticalFunctions()=default;

}
