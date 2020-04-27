#include "mms_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "MMS.H"
#include "MMSForcing.H"
#include "icns/icns.H"
#include "icns/icns_ops.H"

namespace amr_wind_tests {

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(MMSMeshTest, mms_forcing)
{
    constexpr amrex::Real tol = 1.0e-12;

    // Initialize parameters
    utils::populate_mms_params();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    for (auto& pp : sim().physics()) pp->post_init_actions();

    auto fields = ICNSFields(sim().repo())(sim().time(), 0);
    auto& src_term = fields.src_term;
    amr_wind::pde::icns::mms::MMSForcing mmsforcing(sim());

    amrex::Real min_golds[AMREX_SPACEDIM] = {
        -2.1397143441391857, -2.5061563892200622, -2.6756003260809429};
    amrex::Real max_golds[AMREX_SPACEDIM] = {
        2.0381534755116628, 2.2014865191023762, 2.4125363807493985};
    src_term.setVal(0.0);

    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        mmsforcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, min_golds[i], tol);
        EXPECT_NEAR(max_val, max_golds[i], tol);
    }
}
} // namespace amr_wind_tests
