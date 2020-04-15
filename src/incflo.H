#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_NodalProjector.H>

#include <DiffusionTensorOp.H>
#include <DiffusionScalarOp.H>

#include "incflo_enums.H"
#include "SimTime.H"
#include "FieldRepo.H"

namespace amr_wind {
namespace pde {
class PDEBase;
}
class Physics;
class RefinementCriteria;
}

class incflo : public amrex::AmrCore
{
public:

    friend DiffusionTensorOp;
    friend DiffusionScalarOp;

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

    void set_abl_friction_vels(amrex::Real uground, amrex::Real ustar)
    {
        m_velocity_mean_ground = uground;
        m_utau_mean_ground = ustar;
    }

    const amr_wind::SimTime& time() const { return m_time; }
    amr_wind::FieldRepo& repo() { return m_repo; }
    const amr_wind::FieldRepo& repo() const { return m_repo; }
    amr_wind::Field& velocity() const { return m_repo.get_field("velocity"); }
    amr_wind::Field& density() const { return m_repo.get_field("density"); }
    amr_wind::Field& tracer() const { return m_repo.get_field("temperature"); }
    amr_wind::Field& temperature() const { return m_repo.get_field("temperature"); }
    amr_wind::Field& pressure() const { return m_repo.get_field("p"); }
    amr_wind::Field& grad_p() const { return m_repo.get_field("gp"); }

public: // for cuda

    void ComputeDt (bool explicit_diffusion);

    void compute_tra_forces (amrex::Vector<amrex::MultiFab*> const& tra_forces,
                             amrex::Vector<amrex::MultiFab const*> const& density);
    void compute_vel_forces (amrex::Vector<amrex::MultiFab*> const& vel_forces,
                             amrex::Vector<amrex::MultiFab const*> const& velocity,
                             amrex::Vector<amrex::MultiFab const*> const& density,
                             amrex::Vector<amrex::MultiFab const*> const& tracer);
    void compute_vel_forces_on_level ( int lev,
                                       amrex::MultiFab& vel_forces,
                                       const amrex::MultiFab& density,
                                       const amrex::MultiFab& tracer);
    void compute_vel_pressure_terms(int lev, amrex::MultiFab& vel_forces,
                                    const amrex::MultiFab& density);


    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////

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

    bool m_wall_model_flag = false; /** if a wall model is detected on 1 of 6 boundaries it is set to true */

    ///////////////////////////////////////////////////////////////////////////
    //
    // prob
    //
    ///////////////////////////////////////////////////////////////////////////

    void prob_init_fluid (int lev);
    void prob_set_inflow_velocity (int grid_id, amrex::Orientation ori, amrex::Box const& bx,
                                   amrex::Array4<amrex::Real> const& v, int lev, amrex::Real time);

#include "incflo_prob_I.H"

    ///////////////////////////////////////////////////////////////////////////
    //
    // projection
    //
    ///////////////////////////////////////////////////////////////////////////

    void ApplyProjection(amrex::Vector<amrex::MultiFab const*> density,
                         amrex::Real time, amrex::Real scaling_factor, bool incremental);

    ///////////////////////////////////////////////////////////////////////////
    //
    // rheology
    //
    ///////////////////////////////////////////////////////////////////////////

    virtual void compute_viscosity (amrex::Vector<amrex::MultiFab*> const& vel_eta,
                                    amrex::Vector<amrex::MultiFab*> const& tra_eta,
                                    amrex::Vector<amrex::MultiFab const*> const& rho,
                                    amrex::Vector<amrex::MultiFab const*> const& vel,
                                    amrex::Vector<amrex::MultiFab const*> const& tra,
                                    amrex::Real time, int nghost);

    // TODO: Transition
    void declare_fields();

    void init_field_bcs();

private:

    //
    // member variables
    //

    amr_wind::SimTime m_time;
    amr_wind::FieldRepo m_repo;

    std::unique_ptr<amr_wind::pde::PDEBase> m_icns;
    amrex::Vector<std::unique_ptr<amr_wind::pde::PDEBase>> m_scalar_eqns;

    amrex::Vector<std::unique_ptr<amr_wind::Physics>> m_physics;
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

    // AMR / refinement settings
    int m_refine_cutcells = 1;

    amrex::Real m_p000 = 0.0; // background pressure at node (0,0,0).
    amrex::Vector<amrex::Real> m_gp0{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_delp{{0.0, 0.0, 0.0}};

    // Scalar diffusive coefficient
    amrex::Vector<amrex::Real> m_mu_s;

    // Density (if constant)
    amrex::Real m_ro_0 = 1.0;

    bool m_constant_density = true;
    bool m_advect_tracer    = false;
    bool m_test_tracer_conservation = false;

    // Default to MOL (not Godunov)
    bool m_use_godunov         = false;

    // If using Godunov, default to PPM (not PLM)
    bool m_godunov_ppm         = true;

    // If using Godunov, when to add the forcing terms in 
    //    the construction of the "trans" velocities
    bool m_godunov_use_forces_in_trans = false;

    DiffusionType m_diff_type = DiffusionType::Implicit;

    // Fluid properties
    FluidModel m_fluid_model;
    amrex::Real m_mu = 1.0;
    amrex::Real m_Smagorinsky_Lilly_SGS_constant = 0.135;
     
    // fixme this is sloppy and are being used as global variables...
    amrex::Real m_velocity_mean_ground;
    amrex::Real m_utau_mean_ground;
    amrex::Real m_ground_height;

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

    int m_ntrac = 1;

    std::unique_ptr<DiffusionTensorOp> m_diffusion_tensor_op;
    std::unique_ptr<DiffusionScalarOp> m_diffusion_scalar_op;

    //
    // end of member variables
    //

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return m_repo.factory(lev); }

    // Number of ghost cells for field arrays.
    int nghost_state () const {
        return (m_use_godunov) ? 3 : 2;
    }

    int nghost_force () const {
        return (m_use_godunov) ? 1 : 0;
    }

    int nghost_mac () const {
        return (m_use_godunov) ? 1 : 0;
    }

    bool need_divtau () const {
        return !(!m_use_godunov and DiffusionType::Implicit == m_diff_type);
    }

    DiffusionTensorOp* get_diffusion_tensor_op ();
    DiffusionScalarOp* get_diffusion_scalar_op ();

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_projection_bc (amrex::Orientation::Side side) const noexcept;

    void Advance ();
    void AverageDown ();
    void AverageDownTo (int crse_lev);

    ///////////////////////////////////////////////////////////////////////////
    //
    // rheology
    //
    ///////////////////////////////////////////////////////////////////////////

    void ReadRheologyParameters ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // setup
    //
    ///////////////////////////////////////////////////////////////////////////

    void set_background_pressure ();
    void ReadParameters ();
    void ReadIOParameters ();
    void ResizeArrays (); // Resize arrays to fit (up to) max_level + 1 AMR levels
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

    void PrintMaxValues (amrex::Real time);
    void PrintMaxVel (int lev);
    void PrintMaxGp (int lev);
    void CheckForNans (int lev);

};

#endif
