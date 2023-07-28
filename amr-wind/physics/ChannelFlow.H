#ifndef ChannelFlow_H
#define ChannelFlow_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/boundary_conditions/wall_models/WallFunction.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include <AMReX_REAL.H>

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
analytical_smagorinsky_profile(
    const amrex::Real h,
    const amrex::Real Cs,
    const amrex::Real dx,
    const amrex::Real rho,
    const amrex::Real mu,
    const amrex::Real dpdx,
    const amrex::Real C0,
    const amrex::Real C1)
{
    const amrex::Real Cs2 = Cs * Cs;
    const amrex::Real dx2 = dx * dx;
    const amrex::Real Cs2dx2 = Cs2 * dx2;
    return -mu / (rho * Cs2dx2) * h +
           Cs2dx2 / (3 * dpdx) *
               std::pow(
                   (2.0 / Cs2dx2 * dpdx * h +
                    (mu / (rho * Cs2dx2)) * (mu / (rho * Cs2dx2)) + C1),
                   3.0 / 2.0) +
           C0;
}

} // namespace

namespace amr_wind::channel_flow {

/** Channel Flow physics
 *  \ingroup physics
 */
class ChannelFlow : public Physics::Register<ChannelFlow>
{
public:
    static std::string identifier() { return "ChannelFlow"; }

    explicit ChannelFlow(CFDSim& sim);

    ~ChannelFlow() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    //! Use Reichardt function to initialize velocity profile
    template <typename IndexSelector>
    void initialize_fields(
        int level,
        const amrex::Geometry& geom,
        const IndexSelector& idxOp,
        const int n_idx);

    template <typename IndexSelector>
    amrex::Real compute_error(const IndexSelector& idxOp);

    amrex::Real compute_analytical_smagorinsky_error();

    void output_error();

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    //! CFD simulation controller instance
    CFDSim& m_sim;

    const amr_wind::SimTime& m_time;
    FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    WallFunction m_wall_func;

    //! Wall normal directon - Default y direction
    int m_norm_dir{1};

    //! initial density value
    amrex::Real m_rho{1.0};

    //! viscosity
    amrex::Real m_mu{1.0};

    //! Re_tau
    amrex::Real m_re_tau{1000.0};

    //! u_tau
    amrex::Real m_utau{0.01};

    //! y_tau
    amrex::Real m_ytau{0.001};

    //! initial tke value
    amrex::Real m_tke0{0.0};

    //! initial sdr value
    amrex::Real m_sdr0{1000.0};

    //! flag for perturbations of the initial condition
    bool m_perturb_vel{false};

    //! perturbation period (y-direction)
    amrex::Real m_perturb_y_period{1.0};

    //! perturbation period (z-direction)
    amrex::Real m_perturb_z_period{1.0};

    //! perturbation factor (fraction of utau)
    amrex::Real m_perturb_fac{0.1};

    //! Von-Karman constant
    amrex::Real m_kappa{0.41};

    //! Turbulence model
    std::string m_turbulence_model;

    bool m_laminar{false};

    bool m_mesh_mapping{false};

    amrex::Real m_mean_vel;

    bool m_half{false};

    //! flag for analytical smagorinsky test
    bool m_analytical_smagorinsky_test{false};

    //! Analytical Smagorinsky first coefficient
    amrex::Real m_C0{0.0};

    //! Analytical Smagorinsky second coefficient
    amrex::Real m_C1{0.0};

    //! Body forcing (x direction)
    amrex::Real m_dpdx{0.0};

    //! direction of mean velocity
    int m_mean_vel_dir{0};

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"channel_flow.log"};
};
} // namespace amr_wind::channel_flow

#endif /* ChannelFlow_H */
