#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

AnalyticalFunctions::AnalyticalFunctions(
    int n,
    amrex::Box bx)
    : ncells_(n),
      dx_(0.1),
      dy_(0.2 + 0.01*amrex::Random()),
      dz_(0.3),
      scalar_(bx,1),
      scalargrad_(bx,AMREX_SPACEDIM),
      vector_(bx,AMREX_SPACEDIM),
      vectorgrad_(bx,AMREX_SPACEDIM*AMREX_SPACEDIM)
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


LinearAnalyticalFunctions::LinearAnalyticalFunctions(
    int n,
    amrex::Box bx)
    : AnalyticalFunctions(n,bx)
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    auto vec = vector_.array();
    auto gradvec = vectorgrad_.array();
    for(int i = 0; i <= ncells_+1; ++i){
        for(int j = 0; j <= ncells_+1; ++j){
            for(int k = 0; k <= ncells_+1; ++k){
                phi(i,j,k) = 7.8*x_[i]*y_[j]*z_[k] + 3.14;
                gradphi(i,j,k,0) = 7.8 * y_[j] * z_[k];
                gradphi(i,j,k,1) = 7.8 * x_[i] * z_[k];
                gradphi(i,j,k,2) = 7.8 * x_[i] * y_[j];
                vec(i,j,k,0) = 7.8*x_[i]*y_[j]*z_[k] + 3.14;
                vec(i,j,k,1) = 5.1*x_[i]*y_[j]*z_[k] + 9.7123;
                vec(i,j,k,2) = 1.923*x_[i]*y_[j]*z_[k] + 6.19826;
                gradvec(i,j,k,0) = 7.8 * y_[j] * z_[k];
                gradvec(i,j,k,1) = 7.8 * x_[i] * z_[k];
                gradvec(i,j,k,2) = 7.8 * x_[i] * y_[j];
                gradvec(i,j,k,3) = 5.1 * y_[j] * z_[k];
                gradvec(i,j,k,4) = 5.1 * x_[i] * z_[k];
                gradvec(i,j,k,5) = 5.1 * x_[i] * y_[j];
                gradvec(i,j,k,6) = 1.923 * y_[j] * z_[k];
                gradvec(i,j,k,7) = 1.923 * x_[i] * z_[k];
                gradvec(i,j,k,8) = 1.923 * x_[i] * y_[j];
            }
        }
    }
}
LinearAnalyticalFunctions::~LinearAnalyticalFunctions()=default;

QuadraticAnalyticalFunctions::QuadraticAnalyticalFunctions(
    int n,
    amrex::Box bx)
    : AnalyticalFunctions(n,bx)
{
    // Define scalar fields and gradient
    auto phi = scalar_.array();
    auto gradphi = scalargrad_.array();
    auto vec = vector_.array();
    auto gradvec = vectorgrad_.array();
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
                
                vec(i,j,k,0) = 5.0*x*y*z + 4.0*x*x + 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z + 9.4*x*y + 3.4*x+ 2.0*x + 3.14;
                gradvec(i,j,k,0) = 5.0*y*z + 8.0*x + 3.8*z + 9.4*y + 3.4 + 2.0; 
                gradvec(i,j,k,1) = 5.0*x*z + 6.0*y + 1.3*z + 9.4*x; 
                gradvec(i,j,k,2) = 5.0*x*y -4.6*z + 1.3*y + 3.8*x;
                
                vec(i,j,k,1) = 8.265*x*y*z + 1.924*x*x + 0.923*y*y - 8.65*z*z + 2.834*y*z + 9.812*x*z + 4.12*x*y + 1.0921*x + 131.0;
                gradvec(i,j,k,3) = 8.265*y*z + 3.848*x + 9.812*z + 4.12*y + 1.0921; 
                gradvec(i,j,k,4) = 8.265*x*z + 1.846*y + 2.834*z + 4.12*x; 
                gradvec(i,j,k,5) = 8.265*x*y - 17.3*z + 2.834*y + 9.812*x;
                
                vec(i,j,k,2) = 0.2994917334*x*y*z
                    + 0.8228505242180528*x*x + 0.1205690389747357*y*y
                    - 0.9604194947750178*z*z + 0.37039876870122856*y*z
                    + 0.6241061971479255*x*z + 0.7511807179790003*x*y;
                gradvec(i,j,k,6) = 0.2994917334*y*z + 2.0*0.8228505242180528*x
                    + 0.6241061971479255*z + 0.7511807179790003*y; 
                gradvec(i,j,k,7) = 0.2994917334*x*z + 2.0*0.1205690389747357*y
                    + 0.37039876870122856*z + 0.7511807179790003*x; 
                gradvec(i,j,k,8) = 0.2994917334*x*y - 2.0*0.9604194947750178*z
                    + 0.37039876870122856*y + 0.6241061971479255*x;
                
            }
        }
    }
}

QuadraticAnalyticalFunctions::~QuadraticAnalyticalFunctions()=default;



}
