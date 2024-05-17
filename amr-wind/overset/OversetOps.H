#ifndef OVERSETOPS_H
#define OVERSETOPS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind {

class OversetOps
{
public:
    OversetOps() = default;

    void initialize(CFDSim& sim);
    void pre_advance_work();
    void post_advance_work();

    void update_gradp();

private:
    // Functions called within public functions
    void parameter_output();
    void sharpen_nalu_data();
    void set_hydrostatic_gradp();
    void replace_masked_gradp();

    // Check for multiphase sim
    bool m_vof_exists{false};
    // Check for perturbational pressure
    // (will be removed soon)
    bool m_perturb_p{false};

    // This is the only option for now
    const bool m_use_hs_pgrad{true};

    // Coupling options
    bool m_disable_nodal_proj{false};
    bool m_disable_mac_proj{false};
    bool m_replace_gp{false};

    // Verbosity
    int m_verbose{0};

    // Reinitialization parameters
    int m_niterations{10};
    int m_calcconvint{1}; // calctolniter, cconv
    amrex::Real m_tol = 1e-12;
    amrex::Real m_rlscale = 1.5;
    amrex::Real m_margin = 0.1;
    amrex::Real m_target_cutoff = 0.0; // proc_tgvof_tol

    // Tolerance for VOF-related bound checks
    const amrex::Real m_vof_tol = 1e-12;
    // Small number for approximate signed distance function
    const amrex::Real m_asdf_tiny = 1e-12;

    // Pointer for pressure gradient copy field
    amr_wind::Field* m_gp_copy{nullptr};
    // Pointer for MultiPhase physics
    amr_wind::MultiPhase* m_mphase{nullptr};

    CFDSim* m_sim_ptr;
};

} // namespace amr_wind

#endif /* OVERSETOPS_H */