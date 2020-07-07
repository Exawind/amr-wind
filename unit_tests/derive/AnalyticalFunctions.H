#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/derive/derive_K.H"

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
    amrex::FArrayBox scalarlaplace_;
    amrex::FArrayBox vector_;
    amrex::FArrayBox vectorgrad_;
    amrex::FArrayBox vectorlaplace_;
    amrex::FArrayBox curvature_;
    amrex::FArrayBox strainrate_;
};

class LinearAnalyticalFunctions
    : public AnalyticalFunctions
{
public: 
    LinearAnalyticalFunctions(int n, amrex::Box);
    ~LinearAnalyticalFunctions();
};
    
class QuadraticAnalyticalFunctions :
        public AnalyticalFunctions
{
public: 
    QuadraticAnalyticalFunctions(int n, amrex::Box);
    ~QuadraticAnalyticalFunctions();
};

class SphereAnalyticalFunction :
        public AnalyticalFunctions
{
public:
    SphereAnalyticalFunction(int n, amrex::Box);
    ~SphereAnalyticalFunction();
};

}
