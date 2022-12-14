#ifndef ChannelFlow_H
#define ChannelFlow_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include <AMReX_REAL.H>

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

    void output_error();

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

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

    //! Von-Karman constant
    amrex::Real m_kappa{0.41};

    bool m_laminar{false};

    bool m_mesh_mapping{false};

    amrex::Real m_mean_vel;

    //! direction of mean velocity
    int m_mean_vel_dir{0};

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"channel_flow.log"};
};
} // namespace amr_wind::channel_flow

#endif /* ChannelFlow_H */
