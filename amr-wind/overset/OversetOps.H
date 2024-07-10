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
    void parameter_output() const;
    void sharpen_nalu_data();
    void set_hydrostatic_gradp();
    void replace_masked_gradp();

    // Check for multiphase sim
    bool m_vof_exists{false};
    // Check for perturbational pressure
    // (will be removed soon)
    bool m_perturb_p{false};

    // To avoid using sharpened pressure, use hydrostatic pressure
    bool m_use_hydrostatic_gradp{false};

    // Coupling options
    bool m_disable_nodal_proj{false};
    bool m_disable_mac_proj{false};
    bool m_replace_gradp_postsolve{false};

    // Verbosity
    int m_verbose{0};

    // Reinitialization parameters
    int m_n_iterations{10};
    int m_calc_convg_interval{1}; // calctolniter, cconv
    amrex::Real m_convg_tol = 1e-12;
    amrex::Real m_relative_length_scale = 1.5;
    amrex::Real m_upw_margin = 0.1;
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