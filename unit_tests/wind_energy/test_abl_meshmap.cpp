#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/mesh_mapping_models/Geom2ConstantScaling.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind_tests {
namespace {} // namespace

TEST_F(ABLMeshTest, abl_map_initialization)
{
    const amrex::Real tol = 1.0e-12;

    amrex::Real trans_loc = 0.20;
    amrex::Real trans_width = 0.1;
    amrex::Real sratio = 1.1;
    amrex::Real delta0 = 0.1;
    amrex::Real len = 1.0;

    // Test the mapping coord and fac functions
    auto g_mesh_map = amr_wind::geom2constant_map::Geom2ConstantScaling();
    amrex::Real xcoord = g_mesh_map.dump_coord(
        19.0, trans_loc, trans_width, sratio, delta0, len, 1);
    amrex::Real dxcoord = g_mesh_map.dump_fac(
        19.0, trans_loc, trans_width, sratio, delta0, len, 1);
    EXPECT_NEAR(xcoord, 1.9, tol);
    EXPECT_NEAR(dxcoord, 0.1, tol);

    populate_parameters();
    utils::populate_abl_params();

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 48}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0, 120.0, 48.0}};
        pp.addarr("prob_hi", probhi);
    }

    // Initial conditions (Mesh map)
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<std::string> mapping{"Geom2ConstantScaling"};
        pp.addarr("mesh_mapping", mapping);
    }

    // Initial conditions (Scaling)
    {
        amrex::ParmParse pp("Geom2ConstantScaling");
        amrex::Vector<amrex::Real> pp_sratio{{1.1, 1.1, 1.1}};
        amrex::Vector<amrex::Real> pp_delta0{{0.1, 0.1, 1000.0 / 48.0}};
        amrex::Vector<amrex::Real> pp_transwidth{{0.1, 0.1, 0.1}};
        amrex::Vector<amrex::Real> pp_transloc{{0.2, 0.2, 50.0}};
        amrex::Vector<int> pp_do_map{{0, 0, 1}};

        pp.addarr("sratio", pp_sratio);
        pp.addarr("delta0", pp_delta0);
        pp.addarr("transwidth", pp_transwidth);
        pp.addarr("translocation", pp_transloc);
        pp.addarr("do_map", pp_do_map);
    }

    initialize_mesh();

    sim().activate_mesh_map();

    const int nlevels = mesh().num_levels();

    if (sim().has_mesh_mapping()) {
        for (int lev = 0; lev < nlevels; lev++) {
            sim().mesh_mapping()->create_map(lev, mesh().Geom(lev));
        }
    }
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");
    auto& temperaturef = frepo.declare_field("temperature");

    auto velocity = velocityf.vec_ptrs();
    auto density = densityf.vec_ptrs();
    auto temperature = temperaturef.vec_ptrs();

    amr_wind::ABLFieldInit ablinit(sim());
    run_algorithm(
        mesh().num_levels(), density,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            auto rho = density[lev]->array(mfi);
            auto theta = temperature[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            ablinit(lev, mfi, bx, mesh().Geom(lev), vel, rho, theta);
        });

    const int level = 0;
    const int m_direction = 2;
    if (sim().has_mesh_mapping()) {
        const auto& velocity_f = sim().repo().get_field("velocity");
        auto& vel = velocity_f(level);
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell(0.0, AMREX_SPACEDIM);
        pp.queryarr("n_cell", ncell, 0, AMREX_SPACEDIM);
        amrex::Vector<amrex::Real> zpts;
        amrex::Vector<amrex::Real> zfac;
        zpts.resize(ncell[2]);
        zfac.resize(ncell[2]);

        amr_wind::Field const* nu_coord_cc =
            &(sim().repo().get_field("non_uniform_coord_cc"));
        amr_wind::Field const* mesh_scale_fac_cc =
            &(sim().repo().get_field("mesh_scaling_factor_cc"));

        // Loop through and sum over all points on the lower surface
        for (amrex::MFIter mfi(vel); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            amrex::Array4<amrex::Real const> nu_cc =
                ((*nu_coord_cc)(level).array(mfi));
            amrex::Array4<amrex::Real const> fac_cc =
                ((*mesh_scale_fac_cc)(level).array(mfi));

            amrex::Loop(vbx, [=, &zpts, &zfac](int i, int j, int k) noexcept {
                if ((i == 0) && (j == 0)) {
                    zpts[k] = nu_cc(i, j, k, m_direction);
                    zfac[k] = fac_cc(i, j, k, m_direction);
                }
            });
        }
        EXPECT_NEAR(zpts[46], 968.75, tol);
        EXPECT_NEAR(zfac[46], 1000.0 / 48.0, tol);
    } // if (sim().has_mesh_mapping())
}

} // namespace amr_wind_tests
