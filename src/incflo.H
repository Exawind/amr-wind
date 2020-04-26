#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_NodalProjector.H>

#include "incflo_enums.H"
#include "CFDSim.H"
#include "SimTime.H"
#include "FieldRepo.H"

namespace amr_wind {
namespace pde {
class PDEBase;
}
class RefinementCriteria;
}

class incflo : public amrex::AmrCore
{
public:
    incflo ();
    virtual ~incflo ();

    // Initialize multilevel AMR data
    void InitData ();

    // Evolve solution to final time through repeated calls to Advance()
    void Evolve ();

    // Tag cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping
    // Only used during initialization
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time,
                                          const amrex::BoxArray& new_grids,
                                          const amrex::DistributionMapping& new_dmap) override;

    // Make a new level using provided BoxArray and DistributionMapping,
    // and fill with interpolated coarse level data
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping,
    // and fill with existing fine and coarse data
    virtual void RemakeLevel (int lev, amrex::Real time,
                              const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) override;

    // Delete level data
    virtual void ClearLevel (int lev) override;

    amr_wind::CFDSim& sim() { return m_sim; }
    const amr_wind::SimTime& time() const { return m_time; }
    amr_wind::FieldRepo& repo() { return m_repo; }
    const amr_wind::FieldRepo& repo() const { return m_repo; }

    amr_wind::pde::PDEBase& icns() { return m_sim.pde_manager().icns(); }
    amr_wind::pde::PDEMgr::TypeVector& scalar_eqns()
    { return m_sim.pde_manager().scalar_eqns(); }

    amr_wind::Field& velocity() const { return m_repo.get_field("velocity"); }
    amr_wind::Field& density() const { return m_repo.get_field("density"); }
    amr_wind::Field& temperature() const { return m_repo.get_field("temperature"); }
    amr_wind::Field& pressure() const { return m_repo.get_field("p"); }
    amr_wind::Field& grad_p() const { return m_repo.get_field("gp"); }

public: // for cuda

    void ComputeDt (bool explicit_diffusion);

    void set_inflow_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int nghost);

    ///////////////////////////////////////////////////////////////////////////
    //
    // convection
    //
    ///////////////////////////////////////////////////////////////////////////

    void ApplyPredictor(bool incremental_projection = false);
    void ApplyCorrector();

    ///////////////////////////////////////////////////////////////////////////
    //
    // derive
    //
    ///////////////////////////////////////////////////////////////////////////

    void ComputeStrainrate (amrex::Real time);
    void ComputeVorticity (int lev, amrex::Real time, amrex::MultiFab& vort,
                           amrex::MultiFab const& vel);
    void ComputeDivU (amrex::Real time);
    amrex::Real ComputeKineticEnergy () const;

    bool m_wall_model_flag = false; /** if a wall model is detected on at least 1 of 6 boundaries it is set to true */

    ///////////////////////////////////////////////////////////////////////////
    //
    // prob
    //
    ///////////////////////////////////////////////////////////////////////////

    void prob_init_fluid (int lev);
    void prob_set_inflow_velocity (int grid_id, amrex::Orientation ori, amrex::Box const& bx,
                                   amrex::Array4<amrex::Real> const& v, int lev, amrex::Real time,
                                   amrex::Real bcv);

#include "incflo_prob_I.H"

    ///////////////////////////////////////////////////////////////////////////
    //
    // projection
    //
    ///////////////////////////////////////////////////////////////////////////

    void ApplyProjection(amrex::Vector<amrex::MultiFab const*> density,
                         amrex::Real time, amrex::Real scaling_factor, bool incremental);


    // TODO: Transition
    void declare_fields();

    void init_field_bcs();

private:

    //
    // member variables
    //

    amr_wind::CFDSim m_sim;
    amr_wind::SimTime& m_time;
    amr_wind::FieldRepo& m_repo;

    amrex::Vector<std::unique_ptr<amr_wind::RefinementCriteria>> m_refine_criteria;

    // Be verbose?
    int m_verbose = 0;

    // Member variables for initial conditions
    int m_probtype = 0;
    amrex::Real m_ic_u = 0.0;
    amrex::Real m_ic_v = 0.0;
    amrex::Real m_ic_w = 0.0;
    amrex::Real m_ic_p = 0.0;

    // Initial projection / iterations
    bool m_do_initial_proj    = true;
    int  m_initial_iterations = 3;

    amrex::Real m_p000 = 0.0; // background pressure at node (0,0,0).
    amrex::Vector<amrex::Real> m_gp0{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_delp{{0.0, 0.0, 0.0}};

    // Density (if constant)
    amrex::Real m_ro_0 = 1.0;

    bool m_constant_density = true;
    bool m_test_tracer_conservation = false;

    // Default to MOL (not Godunov)
    bool m_use_godunov         = false;

    DiffusionType m_diff_type = DiffusionType::Implicit;

    // Fluid properties
    amrex::Real m_mu = 1.0;

    int m_last_plt = -1;
    std::string m_plot_file{"plt"};

    int m_last_chk = -1;
    int m_KE_int = -1;

    std::string m_check_file{"chk"};
    std::string m_restart_file{""};

    // Flags for saving fluid data in plot files
    int m_plt_velx        = 1;
    int m_plt_vely        = 1;
    int m_plt_velz        = 1;
    int m_plt_gpx         = 1;
    int m_plt_gpy         = 1;
    int m_plt_gpz         = 1;
    int m_plt_rho         = 1;
    int m_plt_tracer      = 1;
    int m_plt_p           = 0;
    int m_plt_eta         = 0;
    int m_plt_vort        = 1;
    int m_plt_forcing     = 0;
    int m_plt_strainrate  = 0;
    int m_plt_stress      = 0;
    int m_plt_divu        = 0;
    int m_plt_vfrac       = 1;

    //
    // end of member variables
    //

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return m_repo.factory(lev); }

    // Number of ghost cells for force field arrays.
    int nghost_force () const {
        return (m_use_godunov) ? 1 : 0;
    }

    bool need_divtau () const {
        return !(!m_use_godunov and DiffusionType::Implicit == m_diff_type);
    }

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_projection_bc (amrex::Orientation::Side side) const noexcept;

    void Advance ();
    void AverageDown ();
    void AverageDownTo (int crse_lev);

    ///////////////////////////////////////////////////////////////////////////
    //
    // setup
    //
    ///////////////////////////////////////////////////////////////////////////

    void set_background_pressure ();
    void ReadParameters ();
    void ReadIOParameters ();
    void InitialProjection ();
    void InitialIterations ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // utilities
    //
    ///////////////////////////////////////////////////////////////////////////

    void WriteHeader (const std::string& name, bool is_checkpoint) const;
    void WriteJobInfo (const std::string& dir) const;
    void WriteCheckPointFile () const;
    void WritePlotFile ();
    void ReadCheckpointFile ();

    void PrintMaxValues (const std::string& header);
    void PrintMaxVel (int lev);
    void PrintMaxGp (int lev);
    void CheckForNans (int lev);

};

#endif
