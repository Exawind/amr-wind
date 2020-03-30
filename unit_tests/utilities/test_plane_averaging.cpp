#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "PlaneAveraging.H"
#include "trig_ops.H"

namespace amr_wind_tests {

class PlaneAveragingTest : public MeshTest
{
public:
    void test_dir(int);
};

TEST_F(PlaneAveragingTest, test_constant)
{
    constexpr double tol = 1.0e-12;
    constexpr amrex::Real u0 = 2.3, v0 = 3.5, w0 = 5.6, t0 = 3.2;

    populate_parameters();
    initialize_mesh();

    auto velocity = mesh().declare_field("velocity", 3);
    auto tracer = mesh().declare_field("tracer");
    
    // initialize level 0 to a constant
    velocity[0]->setVal(u0,0,1);
    velocity[0]->setVal(v0,1,1);
    velocity[0]->setVal(w0,2,1);
    tracer[0]->setVal(t0);
    
    const auto& problo = mesh().Geom(0).ProbLoArray();
    const auto& probhi = mesh().Geom(0).ProbHiArray();
    
    // test the average of a constant is the same constant
    for(int dir=0;dir<3;++dir){
        
        PlaneAveraging pa(mesh().Geom(), velocity, tracer, dir);
        
        amrex::Real z = 0.5*(problo[dir] + probhi[dir]);

        amrex::Real u = pa.line_velocity_xdir(z);
        amrex::Real v = pa.line_velocity_ydir(z);
        amrex::Real w = pa.line_velocity_zdir(z);

        EXPECT_NEAR(u0, u, tol);
        EXPECT_NEAR(v0, v, tol);
        EXPECT_NEAR(w0, w, tol);
    }
    
}
   
TEST_F(PlaneAveragingTest, test_linear)
{
    
        constexpr double tol = 1.0e-12;
        constexpr amrex::Real u0 = 1.0, v0 = 3.5, w0 = 5.6, t0 = 3.2;

        populate_parameters();
        initialize_mesh();

        auto velocity = mesh().declare_field("velocity", 3);
        auto tracer = mesh().declare_field("tracer");
        tracer[0]->setVal(t0);

        constexpr int dir = 2;
    
        const auto& problo = mesh().Geom(0).ProbLoArray();
        const auto& probhi = mesh().Geom(0).ProbHiArray();
        
        
        run_algorithm(mesh().num_levels(), tracer,
            [&](const int lev, const amrex::MFIter& mfi) {
            
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();
            const auto& dx = mesh().Geom(lev).CellSizeArray();
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                
                amrex::Real x[3];
                
                x[0] = problo[0] + (i + 0.5) * dx[0];
                x[1] = problo[1] + (j + 0.5) * dx[1];
                x[2] = problo[2] + (k + 0.5) * dx[2];

                vel(i,j,k,0) = u0*x[dir];
                vel(i,j,k,1) = v0*x[dir];
                vel(i,j,k,2) = w0*x[dir];

                
            });
            
        });
        

        PlaneAveraging pa(mesh().Geom(), velocity, tracer, dir);
        
        constexpr int n = 20;
        const amrex::Real L = probhi[dir] - problo[dir];
        const amrex::Real dx = L/((amrex::Real) n);
        const amrex::Real hchLo = problo[dir] + 0.5*mesh().Geom(0).CellSizeArray()[dir];
        const amrex::Real hchHi = probhi[dir] - 0.5*mesh().Geom(0).CellSizeArray()[dir];

        for(int i=0;i<n;++i){
            const amrex::Real x = problo[dir] + i*dx;
            
            const amrex::Real u = pa.line_velocity_xdir(x);
            const amrex::Real v = pa.line_velocity_ydir(x);
            const amrex::Real w = pa.line_velocity_zdir(x);
            
            if(x < hchLo){
                // test near bottom boundary where solution is not linearly interpolated
                // but instead copied from first cell
                EXPECT_NEAR(u0*hchLo, u, tol);
                EXPECT_NEAR(v0*hchLo, v, tol);
                EXPECT_NEAR(w0*hchLo, w, tol);
            } else if(x > hchHi){
                // test near top boundary where solution is not linearly interpolated
                // but instead copied from last cell
                EXPECT_NEAR(u0*hchHi, u, tol);
                EXPECT_NEAR(v0*hchHi, v, tol);
                EXPECT_NEAR(w0*hchHi, w, tol);
            } else {
                // test linear interpolation
                EXPECT_NEAR(u0*x, u, tol);
                EXPECT_NEAR(v0*x, v, tol);
                EXPECT_NEAR(w0*x, w, tol);
            }
        }
    
}


void PlaneAveragingTest::test_dir(int dir)
{
    
        constexpr double tol = 1.0e-12;
        constexpr amrex::Real u0 = 2.3, v0 = 3.5, w0 = 5.6, t0 = 3.2;

        populate_parameters();
        initialize_mesh();

        auto velocity = mesh().declare_field("velocity", 3);
        auto tracer = mesh().declare_field("tracer");
        tracer[0]->setVal(t0);

        constexpr int periods = 3;
        
        const auto& problo = mesh().Geom(0).ProbLoArray();
        const auto& probhi = mesh().Geom(0).ProbHiArray();
        
        amrex::Real cx[3];
        cx[0] = periods * amr_wind::utils::two_pi() / (probhi[0] - problo[0]);
        cx[1] = periods * amr_wind::utils::two_pi() / (probhi[1] - problo[1]);
        cx[2] = periods * amr_wind::utils::two_pi() / (probhi[2] - problo[2]);
        
        run_algorithm(mesh().num_levels(), tracer,
            [&](const int lev, const amrex::MFIter& mfi) {
            
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();
            const auto& dx = mesh().Geom(lev).CellSizeArray();
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                
                amrex::Real x[3];
                
                x[0] = problo[0] + (i + 0.5) * dx[0];
                x[1] = problo[1] + (j + 0.5) * dx[1];
                x[2] = problo[2] + (k + 0.5) * dx[2];
                
                vel(i,j,k,0) = u0;
                vel(i,j,k,1) = v0;
                vel(i,j,k,2) = w0;

                for(int d=0;d<3;++d){
                    if(d!=dir){
                        vel(i,j,k,0) += std::cos(cx[d]*x[d]);
                        vel(i,j,k,1) += std::sin(cx[d]*x[d]);
                        vel(i,j,k,2) += std::sin(cx[d]*x[d])*cos(cx[d]*x[d]);
                    }
                }
                
            });
            
        });
        

        PlaneAveraging pa(mesh().Geom(), velocity, tracer, dir);
        
        amrex::Real x = 0.5*(problo[dir] + probhi[dir]);
        amrex::Real u = pa.line_velocity_xdir(x);
        amrex::Real v = pa.line_velocity_ydir(x);
        amrex::Real w = pa.line_velocity_zdir(x);

        // test that a periodic function orthogonal to dir
        // averages out to a constant
        EXPECT_NEAR(u0, u, tol);
        EXPECT_NEAR(v0, v, tol);
        EXPECT_NEAR(w0, w, tol);
        
}
    
TEST_F(PlaneAveragingTest, test_xdir)
{
    test_dir(0);
}
TEST_F(PlaneAveragingTest, test_ydir)
{
    test_dir(1);
}
TEST_F(PlaneAveragingTest, test_zdir)
{
    test_dir(2);
}

TEST_F(PlaneAveragingTest, test_leveldata_zdir)
{
    
    constexpr double tol = 1.0e-12;
    constexpr amrex::Real u0 = 2.3, v0 = 3.5, w0 = 5.6, t0 = 3.2;
    constexpr int periods = 3;
    constexpr int dir = 2;
    
    populate_parameters();
    initialize_mesh();

    auto tracer = mesh().declare_field("tracer");
    
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > new_fact(new amrex::FArrayBoxFactory());
    std::unique_ptr<LevelData> ld (new LevelData(tracer[0]->boxArray(), tracer[0]->DistributionMap(), *new_fact, 1, 1,0,1,1));
  
    ld->tracer.setVal(t0);
 
    const auto& problo = mesh().Geom(0).ProbLoArray();
    const auto& probhi = mesh().Geom(0).ProbHiArray();
    
    amrex::Real cx[3];
    cx[0] = periods * amr_wind::utils::two_pi() / (probhi[0] - problo[0]);
    cx[1] = periods * amr_wind::utils::two_pi() / (probhi[1] - problo[1]);
    cx[2] = periods * amr_wind::utils::two_pi() / (probhi[2] - problo[2]);
   
    amrex::Vector<amrex::MultiFab*> r(1);
    r[0] = &(ld->tracer);
    
    run_algorithm(1, r,
        [&](const int lev, const amrex::MFIter& mfi) {
        
        auto vel = ld->velocity.array(mfi);
        const auto& bx = mfi.validbox();
        const auto& dx = mesh().Geom(lev).CellSizeArray();
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            
            amrex::Real x[3];
            
            x[0] = problo[0] + (i + 0.5) * dx[0];
            x[1] = problo[1] + (j + 0.5) * dx[1];
            x[2] = problo[2] + (k + 0.5) * dx[2];
            
            vel(i,j,k,0) = u0;
            vel(i,j,k,1) = v0;
            vel(i,j,k,2) = w0;

            for(int d=0;d<3;++d){
                if(d!=dir){
                    vel(i,j,k,0) += std::cos(cx[d]*x[d]);
                    vel(i,j,k,1) += std::sin(cx[d]*x[d]);
                    vel(i,j,k,2) += std::sin(cx[d]*x[d])*cos(cx[d]*x[d]);
                }
            }
            
        });
        
    });
    
    PlaneAveraging pa(mesh().Geom(0), *ld, dir);
    
    amrex::Real x = 0.5*(problo[dir] + probhi[dir]);
    amrex::Real u = pa.line_velocity_xdir(x);
    amrex::Real v = pa.line_velocity_ydir(x);
    amrex::Real w = pa.line_velocity_zdir(x);

    // test that a periodic function orthogonal to dir
    // averages out to a constant
    EXPECT_NEAR(u0, u, tol);
    EXPECT_NEAR(v0, v, tol);
    EXPECT_NEAR(w0, w, tol);
        
}

}  // amr_wind_tests
