#include "aw_test_utils/MeshTest.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind_tests {

class MultiphaseTest : public MeshTest
{
public:
    void declare_default_fields()
    {
        auto& field_repo = mesh().field_repo();
        auto& vel = field_repo.declare_field("velocity", 3, 1, 2);
        auto& mask_cell = field_repo.declare_int_field("mask_cell", 1, 1);
        auto& umac = field_repo.declare_field(
            "u_mac", 1, 0, 1, amr_wind::FieldLoc::XFACE);
        auto& vmac = field_repo.declare_field(
            "v_mac", 1, 0, 1, amr_wind::FieldLoc::YFACE);
        auto& wmac = field_repo.declare_field(
            "w_mac", 1, 0, 1, amr_wind::FieldLoc::ZFACE);
    }

protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
        amrex::ParmParse pp("incflo");
        pp.add("use_godunov", 1);
        pp.add("use_limiter",0);
        }
        
        {
            amrex::ParmParse pp("time");
            pp.add("stop_time",3.);
            pp.add("max_step", 3000);
            pp.add("fixed_dt", 0.001);
            pp.add("plot_interval",50);
        }
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{64, 64, 64}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 64);
            pp.add("blocking_factor", 64);
            pp.addarr("n_cell", ncell);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{1.0, 1.0, 1.0}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
    amrex::Real TT=3.;

};

TEST_F(MultiphaseTest, levelset_test2)
{
    populate_parameters();
    create_mesh_instance();
    declare_default_fields();
    initialize_mesh();

    // Register the levelset pde
    auto& pde_mgr = mesh().sim().pde_manager();
    pde_mgr.register_icns();        
    auto& levelset_pde=pde_mgr.register_transport_pde("Levelset");

    // Initialize the input/output manager
    auto& iomgr = sim().io_manager();
    iomgr.register_output_var("levelset");
    iomgr.register_output_var("u_mac");
    iomgr.register_output_var("v_mac");
    iomgr.register_output_var("w_mac");

    sim().io_manager().initialize_io();
    auto& time = sim().time();

    auto& field_repo = mesh().field_repo();
    const int nlevels = field_repo.num_active_levels();

    //Declare the initial levelset field
    auto& levelset = field_repo.get_field("levelset").state(amr_wind::FieldState::New);
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& dx = mesh().Geom(lev).CellSizeArray();  
            const auto& problo = mesh().Geom(lev).ProbLoArray();
            const auto& probhi = mesh().Geom(lev).ProbHiArray();
      
            auto phi = levelset(lev).array(mfi);

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i+0.5)*dx[0];
                const amrex::Real y = problo[1] + (j+0.5)*dx[1];
                const amrex::Real z = problo[2] + (k+0.5)*dx[2];
                const amrex::Real x0 = 0.5 * (problo[0] + probhi[0]);
                const amrex::Real y0 = 0.5 * (problo[1] + probhi[1]);
                const amrex::Real z0 = 0.5 * (problo[2] + probhi[2]);
                const amrex::Real R=0.15;
                phi(i, j, k) = R - std::sqrt((x-x0+R) * (x - x0+R) + 
                                   (y - y0+R) * (y - y0+R) + (z - z0+R) * (z - z0+R));
            });
        }  
    }      
    
    
    // Declare mac velocities
    auto& umac = field_repo.declare_field(
        "u_mac", 1, 0, 1, amr_wind::FieldLoc::XFACE);
    auto& vmac = field_repo.declare_field(
        "v_mac", 1, 0, 1, amr_wind::FieldLoc::YFACE);
    auto& wmac = field_repo.declare_field(
        "w_mac", 1, 0, 1, amr_wind::FieldLoc::ZFACE);

    iomgr.write_plot_file();
    while (time.new_timestep()) {

    amrex::Print()<< "Current time is " <<time.current_time()<<std::endl;
        // Assign umac, vmac and wmac velocities
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(umac(lev)); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                const auto& dx = mesh().Geom(lev).CellSizeArray();  
                const auto& problo = mesh().Geom(lev).ProbLoArray();
               
                auto ux = umac(lev).array(mfi);
                auto uy = vmac(lev).array(mfi);
                auto uz = wmac(lev).array(mfi);
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real xface = problo[0] + i*dx[0];
                    const amrex::Real yface = problo[1] + j*dx[1];
                    const amrex::Real zface = problo[2] + k*dx[2];
                    const amrex::Real xc = problo[0] + (i+0.5)*dx[0];
                    const amrex::Real yc = problo[1] + (j+0.5)*dx[1];
                    const amrex::Real zc = problo[2] + (k+0.5)*dx[2];
                    ux(i,j,k) = 2.0*std::sin(M_PI*xface)*std::sin(M_PI*xface)
                                   *std::sin(2.0*M_PI*yc)*std::sin(2.0*M_PI*zc)*std::cos(M_PI*time.current_time()/TT);
                    uy(i,j,k) = -std::sin(M_PI*yface)*std::sin(M_PI*yface)
                                   *std::sin(2.0*M_PI*xc)*std::sin(2.0*M_PI*zc)*std::cos(M_PI*time.current_time()/TT);
                    uz(i,j,k) = -std::sin(M_PI*zface)*std::sin(M_PI*zface)
                                   *std::sin(2.0*M_PI*xc)*std::sin(2.0*M_PI*yc)*std::cos(M_PI*time.current_time()/TT);
                });
            }
        }    
        // Pre-advance 
        levelset_pde.fields().field.advance_states();
        levelset_pde.fields().field.state(amr_wind::FieldState::Old).fillpatch(time.current_time());

        // Do advection
        levelset_pde.compute_advection_term(amr_wind::FieldState::Old);
        levelset_pde.compute_predictor_rhs(DiffusionType::Explicit);
        
        // output plots
        if(time.write_plot_file()){
            iomgr.write_plot_file();
        }
    }
}

} // namespace amr_wind_tests
