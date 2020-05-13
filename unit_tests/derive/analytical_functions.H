#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"

namespace amr_wind_tests {

class AnalyticalFunctions 
{
public:
    AnalyticalFunctions(int n, amrex::Box);

    ~AnalyticalFunctions();

    int ncells_;
    amrex::Real dx_,dy_,dz_;
    amrex::Vector<amrex::Real> x_, y_, z_;
    amrex::FArrayBox scalar_;
    amrex::FArrayBox scalargrad_;

};

AnalyticalFunctions::AnalyticalFunctions(
    int n,
    amrex::Box bx)
  : ncells_(n),
    dx_(0.1),
    dy_(0.2 + 0.01*amrex::Random()),
    dz_(0.3),
    scalar_(bx,1),
    scalargrad_(bx,AMREX_SPACEDIM)
{
    // Coordinates - First and last are on the boundary faces
    x_.resize(ncells_+2);
    y_.resize(ncells_+2);
    z_.resize(ncells_+2);
    x_[0] = 0.3; x_[ncells_+1] = x_[0] + ncells_*dx_;
    y_[0] = 0.2; y_[ncells_+1] = y_[0] + ncells_*dy_;
    z_[0] = 0.1; z_[ncells_+1] = z_[0] + ncells_*dz_;
    for (int i=1; i <= ncells_; ++i) {
        x_[i] = x_[0] + i*dx_ - 0.5*dx_;
        y_[i] = y_[0] + i*dy_ - 0.5*dy_;
        z_[i] = z_[0] + i*dz_ - 0.5*dz_;        
    }
}

AnalyticalFunctions::~AnalyticalFunctions()=default;

class LinearAnalyticalFunctions
    : public AnalyticalFunctions
{
public: 
    LinearAnalyticalFunctions(int n, amrex::Box);
    ~LinearAnalyticalFunctions();
};

LinearAnalyticalFunctions::LinearAnalyticalFunctions(
    int n,
    amrex::Box bx)
    : AnalyticalFunctions(n,bx)
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    for(int i = 0; i <= ncells_+1; ++i){
        for(int j = 0; j <= ncells_+1; ++j){
            for(int k = 0; k <= ncells_+1; ++k){
                phi(i,j,k) = 7.8*x_[i]*y_[j]*z_[k] + 3.14;
                gradphi(i,j,k,0) = 7.8 * y_[j] * z_[k];
                gradphi(i,j,k,1) = 7.8 * x_[i] * z_[k];
                gradphi(i,j,k,2) = 7.8 * x_[i] * y_[j];
            }
        }
    }
}
LinearAnalyticalFunctions::~LinearAnalyticalFunctions()=default;
    
class QuadraticAnalyticalFunctions :
        public AnalyticalFunctions
{
public: 
    QuadraticAnalyticalFunctions(int n, amrex::Box);
    ~QuadraticAnalyticalFunctions();
};

QuadraticAnalyticalFunctions::QuadraticAnalyticalFunctions(
    int n,
    amrex::Box bx)
    : AnalyticalFunctions(n,bx)
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    for(int i = 0; i <= ncells_+1; ++i){
        for(int j = 0; j <= ncells_+1; ++j){
            for(int k = 0; k <= ncells_+1; ++k){
                amrex::Real x = x_[i];
                amrex::Real y = y_[j];
                amrex::Real z = z_[k];
                phi(i,j,k) = 5.0*x*y*z + 4.0*x*x + 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z + 9.4*x*y + 3.4*x+ 2.0*x + 3.14;
                gradphi(i,j,k,0) = 5.0*y*z + 8.0*x + 3.8*z + 9.4*y + 3.4 + 2.0; 
                gradphi(i,j,k,1) = 5.0*x*z + 6.0*y + 1.3*z + 9.4*x; 
                gradphi(i,j,k,2) = 5.0*x*y -4.6*z + 1.3*y + 3.8*x;
            }
        }
    }
}

QuadraticAnalyticalFunctions::~QuadraticAnalyticalFunctions()=default;


}
