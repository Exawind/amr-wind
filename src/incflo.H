#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_NodalProjector.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <DiffusionTensorOp.H>
#include <DiffusionScalarOp.H>

class incflo : public amrex::AmrCore
{
public:

    friend DiffusionTensorOp;
    friend DiffusionScalarOp;

    enum struct FluidModel {
        Newtonian, powerlaw, Bingham, HerschelBulkley, deSouzaMendesDutra
    };

    incflo ();
    virtual ~incflo ();

    // Initialize multilevel AMR data
    void InitData ();

#ifdef AMREX_USE_EB
    void WriteMyEBSurface ();
#endif

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

public: // for cuda

    void ComputeDt (int initialisation, bool explicit_diffusion);
    void compute_forces (amrex::Vector<amrex::MultiFab*> const& vel_forces,
                         amrex::Vector<amrex::MultiFab*> const& tra_forces,
                         amrex::Vector<amrex::MultiFab const*> const& density,
                         amrex::Vector<amrex::MultiFab const*> const& tracer);

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
    void compute_convective_term (amrex::Vector<amrex::MultiFab*> const& conv_u,
                                  amrex::Vector<amrex::MultiFab*> const& conv_r,
                                  amrex::Vector<amrex::MultiFab*> const& conv_t,
                                  amrex::Vector<amrex::MultiFab const*> const& vel,
                                  amrex::Vector<amrex::MultiFab const*> const& density,
                                  amrex::Vector<amrex::MultiFab const*> const& tracer,
//                                amrex::Vector<amrex::MultiFab const*> const& u_mac,
//                                amrex::Vector<amrex::MultiFab const*> const& v_mac,
//                                amrex::Vector<amrex::MultiFab const*> const& w_mac,
                                  amrex::Vector<amrex::MultiFab*> const& u_mac,
                                  amrex::Vector<amrex::MultiFab*> const& v_mac,
                                  amrex::Vector<amrex::MultiFab*> const& w_mac,
                                  amrex::Vector<amrex::MultiFab const*> const& vel_forces,
                                  amrex::Vector<amrex::MultiFab const*> const& tra_forces,
                                  amrex::Real time);

    void compute_convective_term (amrex::Box const& bx, int lev, amrex::MFIter const& mfi,
                                  amrex::Array4<amrex::Real> const& dvdt, // velocity
                                  amrex::Array4<amrex::Real> const& drdt, // density
                                  amrex::Array4<amrex::Real> const& dtdt, // tracer
                                  amrex::Array4<amrex::Real const> const& vel,
                                  amrex::Array4<amrex::Real const> const& rho,
                                  amrex::Array4<amrex::Real const> const& tra,
                                  amrex::Array4<amrex::Real const> const& umac,
                                  amrex::Array4<amrex::Real const> const& vmac,
                                  amrex::Array4<amrex::Real const> const& wmac,
                                  amrex::Array4<amrex::Real const> const& fv,
                                  amrex::Array4<amrex::Real const> const& ft);

    //
    // Godunov
    //
    void predict_godunov (int lev, amrex::Real time, amrex::MultiFab& u_mac,
                          amrex::MultiFab& v_mac, amrex::MultiFab& w_mac,
                          amrex::MultiFab const& vel, amrex::MultiFab const& rho,
                          amrex::MultiFab const& vel_forces);

    void make_ppm_integrals (int lev, amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& Imx,
                             amrex::Array4<amrex::Real> const& Ipx,
                             amrex::Array4<amrex::Real> const& Imy,
                             amrex::Array4<amrex::Real> const& Ipy,
                             amrex::Array4<amrex::Real> const& Imz,
                             amrex::Array4<amrex::Real> const& Ipz,
                             amrex::Array4<amrex::Real const> const& q,
                             amrex::Array4<amrex::Real const> const& vel);

    void make_trans_velocities (int lev, amrex::Box const& xbx,
                                amrex::Box const& ybx, amrex::Box const& zbx,
                                amrex::Array4<amrex::Real> const& u_ad,
                                amrex::Array4<amrex::Real> const& v_ad,
                                amrex::Array4<amrex::Real> const& w_ad,
                                amrex::Array4<amrex::Real const> const& Imx,
                                amrex::Array4<amrex::Real const> const& Ipx,
                                amrex::Array4<amrex::Real const> const& Imy,
                                amrex::Array4<amrex::Real const> const& Ipy,
                                amrex::Array4<amrex::Real const> const& Imz,
                                amrex::Array4<amrex::Real const> const& Ipz,
                                amrex::Array4<amrex::Real const> const& vel,
                                amrex::Array4<amrex::Real const> const& f);

    void predict_godunov (int lev, amrex::Box const& bx, int ncomp,
                          amrex::Box const& xbx, amrex::Box const& ybx, amrex::Box const& zbx,
                          amrex::Array4<amrex::Real> const& qx,
                          amrex::Array4<amrex::Real> const& qy,
                          amrex::Array4<amrex::Real> const& qz,
                          amrex::Array4<amrex::Real const> const& q,
                          amrex::Array4<amrex::Real const> const& u_ad,
                          amrex::Array4<amrex::Real const> const& v_ad,
                          amrex::Array4<amrex::Real const> const& w_ad,
                          amrex::Array4<amrex::Real> const& Imx,
                          amrex::Array4<amrex::Real> const& Ipx,
                          amrex::Array4<amrex::Real> const& Imy,
                          amrex::Array4<amrex::Real> const& Ipy,
                          amrex::Array4<amrex::Real> const& Imz,
                          amrex::Array4<amrex::Real> const& Ipz,
                          amrex::Array4<amrex::Real const> const& f,
                          amrex::Real* p);

    void compute_godunov_advection (int lev, amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& dqdt,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& umac,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::Array4<amrex::Real const> const& wmac,
                                    amrex::Array4<amrex::Real const> const& fq,
                                    amrex::BCRec const* d_bcrec, int const* iconserv,
                                    amrex::Real* p);

    //
    // MOL
    //
    void predict_vels_on_faces (int lev, amrex::Real time, amrex::MultiFab& u_mac,
                                amrex::MultiFab& v_mac, amrex::MultiFab& w_mac,
                                amrex::MultiFab const& vel);

    void predict_vels_on_faces (int lev,
                                amrex::Box const& ubx,
                                amrex::Box const& vbx,
                                amrex::Box const& wbx,
                                amrex::Array4<amrex::Real> const& u,
                                amrex::Array4<amrex::Real> const& v,
                                amrex::Array4<amrex::Real> const& w,
                                amrex::Array4<amrex::Real const> const& vcc);

#ifdef AMREX_USE_EB
    void predict_vels_on_faces_eb (int lev,
                                   amrex::Box const& ccbx,
                                   amrex::Box const& ubx,
                                   amrex::Box const& vbx,
                                   amrex::Box const& wbx,
                                   amrex::Array4<amrex::Real> const& u,
                                   amrex::Array4<amrex::Real> const& v,
                                   amrex::Array4<amrex::Real> const& w,
                                   amrex::Array4<amrex::Real const> const& vcc,
                                   amrex::Array4<amrex::EBCellFlag const> const& flag,
                                   amrex::Array4<amrex::Real const> const& fcx,
                                   amrex::Array4<amrex::Real const> const& fcy,
                                   amrex::Array4<amrex::Real const> const& fcz, 
                                   amrex::Array4<amrex::Real const> const& ccc);
#endif
     
    void incflo_correct_small_cells (amrex::Vector<amrex::MultiFab*      > const& cc_vel,
                                     amrex::Vector<amrex::MultiFab const*> const& u_mac,
                                     amrex::Vector<amrex::MultiFab const*> const& v_mac,
                                     amrex::Vector<amrex::MultiFab const*> const& w_mac);

    void compute_convective_fluxes (int lev, amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& fx,
                                    amrex::Array4<amrex::Real> const& fy,
                                    amrex::Array4<amrex::Real> const& fz,
                                    amrex::Array4<amrex::Real const> const& q,
                                    amrex::Array4<amrex::Real const> const& umac,
                                    amrex::Array4<amrex::Real const> const& vmac,
                                    amrex::Array4<amrex::Real const> const& wmac,
                                    amrex::BCRec const* h_bcrec,
                                    amrex::BCRec const* d_bcrec);

    void compute_convective_rate (int lev, amrex::Box const& bx, int ncomp,
                                  amrex::Array4<amrex::Real> const& dUdt,
                                  amrex::Array4<amrex::Real const> const& fx,
                                  amrex::Array4<amrex::Real const> const& fy,
                                  amrex::Array4<amrex::Real const> const& fz);

#ifdef AMREX_USE_EB
    void compute_convective_fluxes_eb (int lev, amrex::Box const& bx, int ncomp,
                                       amrex::Array4<amrex::Real> const& fx,
                                       amrex::Array4<amrex::Real> const& fy,
                                       amrex::Array4<amrex::Real> const& fz,
                                       amrex::Array4<amrex::Real const> const& q,
                                       amrex::Array4<amrex::Real const> const& umac,
                                       amrex::Array4<amrex::Real const> const& vmac,
                                       amrex::Array4<amrex::Real const> const& wmac,
                                       amrex::BCRec const* h_bcrec,
                                       amrex::BCRec const* d_bcrec,
                                       amrex::Array4<amrex::EBCellFlag const> const& flag,
                                       amrex::Array4<amrex::Real const> const& fcx,
                                       amrex::Array4<amrex::Real const> const& fcy,
                                       amrex::Array4<amrex::Real const> const& fcz,
                                       amrex::Array4<amrex::Real const> const& ccc);

    void compute_convective_rate_eb (int lev, amrex::Box const& bx, int ncomp,
                                     amrex::Array4<amrex::Real> const& dUdt,
                                     amrex::Array4<amrex::Real const> const& fx,
                                     amrex::Array4<amrex::Real const> const& fy,
                                     amrex::Array4<amrex::Real const> const& fz,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     amrex::Array4<amrex::Real const> const& apx,
                                     amrex::Array4<amrex::Real const> const& apy,
                                     amrex::Array4<amrex::Real const> const& apz);

    void redistribute_eb (int lev, amrex::Box const& bx, int ncomp,
                          amrex::Array4<amrex::Real> const& dUdt,
                          amrex::Array4<amrex::Real const> const& dUdt_in,
                          amrex::Array4<amrex::Real> const& scratch,
                          amrex::Array4<amrex::EBCellFlag const> const& flag,
                          amrex::Array4<amrex::Real const> const& vfrac);
#endif

    ///////////////////////////////////////////////////////////////////////////
    //
    // derive
    //
    ///////////////////////////////////////////////////////////////////////////

    void ComputeStrainrate (amrex::Real time);
    void ComputeVorticity (int lev, amrex::Real time, amrex::MultiFab& vort,
                           amrex::MultiFab const& vel);
    void ComputeDivU (amrex::Real time);
    void ComputeDrag ();
    amrex::Real ComputeKineticEnergy () const;

    ///////////////////////////////////////////////////////////////////////////
    //
    // diffusion
    //
    ///////////////////////////////////////////////////////////////////////////

    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
    average_velocity_eta_to_faces (int lev, amrex::MultiFab const& cc_eta) const;

    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
    average_tracer_eta_to_faces (int lev, int comp, amrex::MultiFab const& cc_eta) const;

    void fixup_eta_on_domain_faces (int lev,
                                    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fc,
                                    amrex::MultiFab const& cc) const;

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

    void ApplyProjection(amrex::Real time, amrex::Real scaling_factor, bool incremental);

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

    ///////////////////////////////////////////////////////////////////////////
    //
    // utilities
    //
    ///////////////////////////////////////////////////////////////////////////

    amrex::Real volWgtSum (int lev, const amrex::MultiFab& mf, int comp, bool local=false);
    bool SteadyStateReached ();

private:

    //
    // member variables
    //

    // Be verbose?
    int m_verbose = 0;

    // Member variables for initial conditions
    int m_probtype = 0;
    amrex::Real m_ic_u = 0.0;
    amrex::Real m_ic_v = 0.0;
    amrex::Real m_ic_w = 0.0;
    amrex::Real m_ic_p = 0.0;

    amrex::Vector<amrex::Real> m_t_old;
    amrex::Vector<amrex::Real> m_t_new;

    // Times
    amrex::Real m_cur_time     = 0.0;
    amrex::Real m_dt           = -1.0;
    amrex::Real m_prev_dt      = -1.0;
    amrex::Real m_prev_prev_dt = -1.0;

    // Time step counter
    int m_nstep         = -1;

    // Stop simulation if cur_time reaches stop_time OR nstep reaches max_step
    // OR steady_state = true AND steady_state_tol is reached
    amrex::Real m_stop_time = -1.0;
    int m_max_step = -1;
    bool m_steady_state = false;
    amrex::Real m_steady_state_tol = 1.0e-5;

    // Options to control time stepping
    amrex::Real m_cfl = 0.5;
    amrex::Real m_fixed_dt = -1.;
    amrex::Real m_init_shrink = 0.1;

    // Initial projection / iterations
    bool m_do_initial_proj    = true;
    int  m_initial_iterations = 3;

    // Use Boussinesq approximation for buoyancy?
    bool m_use_boussinesq = false;

    // AMR / refinement settings
    int m_refine_cutcells = 1;
    int m_regrid_int = -1;

    int m_mac_mg_verbose = 0;
    int m_mac_mg_cg_verbose = 0;
    int m_mac_mg_maxiter = 200;
    int m_mac_mg_cg_maxiter = 200;
    int m_mac_mg_max_coarsening_level = 100;
    amrex::Real m_mac_mg_rtol = 1.0e-11;
    amrex::Real m_mac_mg_atol = 1.0e-14;
    std::string m_mac_bottom_solver_type;

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

    bool m_use_godunov         = false;
    bool m_use_forces_in_trans = false;

    enum struct DiffusionType {
        Invalid, Explicit, Crank_Nicolson, Implicit
    };
    DiffusionType m_diff_type = DiffusionType::Implicit;

    // Fluid properties
    FluidModel m_fluid_model;
    amrex::Real m_mu = 1.0;
    amrex::Real m_n_0 = 0.0;
    amrex::Real m_tau_0 = 0.0;
    amrex::Real m_papa_reg = 0.0;
    amrex::Real m_eta_0 = 0.0;

    int m_plot_int = -1;

    // Dump plotfiles at as close as possible to the designated period *without* changing dt
    amrex::Real m_plot_per_approx = -1.0;
    
    // Dump plotfiles at exactcly the designated period by changing dt
    amrex::Real m_plot_per_exact  = -1.0;

    int m_last_plt = -1;
    std::string m_plot_file{"plt"};

    int m_check_int = -1;
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
    int m_plt_strainrate  = 0;
    int m_plt_stress      = 0;
    int m_plt_divu        = 0;
    int m_plt_vfrac       = 1;

    struct LevelData {
        LevelData () = default;
        LevelData (amrex::BoxArray const& ba,
                   amrex::DistributionMapping const& dm,
                   amrex::FabFactory<amrex::FArrayBox> const& fact,
                   int ntrac, int ng_state,
                   bool use_godunov, bool implicit_diffusion, bool advect_tracer);
        // cell-centered multifabs
        amrex::MultiFab velocity;
        amrex::MultiFab velocity_o;
        amrex::MultiFab density;
        amrex::MultiFab density_o;
        amrex::MultiFab tracer;
        amrex::MultiFab tracer_o;
        amrex::MultiFab gp; // pressure gradient
        // nodal multifabs
        amrex::MultiFab p;
        //
        amrex::MultiFab conv_velocity;
        amrex::MultiFab conv_velocity_o;
        amrex::MultiFab conv_density;
        amrex::MultiFab conv_density_o;
        amrex::MultiFab conv_tracer;
        amrex::MultiFab conv_tracer_o;
        //
        amrex::MultiFab divtau;
        amrex::MultiFab divtau_o;
        amrex::MultiFab laps;
        amrex::MultiFab laps_o;
    };

    amrex::Vector<std::unique_ptr<LevelData> > m_leveldata;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

    enum struct BC {
        pressure_inflow, pressure_outflow, mass_inflow, no_slip_wall, slip_wall,
        periodic, undefined
    };

    amrex::GpuArray<BC                         , AMREX_SPACEDIM*2> m_bc_type;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_pressure;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_density;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                               , AMREX_SPACEDIM*2> m_bc_velocity;
    amrex::GpuArray<amrex::Vector<amrex::Real> , AMREX_SPACEDIM*2> m_bc_tracer;

    // amrex::Vector cannot be used on gpu, so ...
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> m_bc_tracer_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_bc_tracer_raii;
    //
    amrex::Vector<amrex::BCRec> m_bcrec_velocity;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_velocity_d;
    amrex::Vector<amrex::BCRec> m_bcrec_density;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_density_d;
    amrex::Vector<amrex::BCRec> m_bcrec_tracer;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_tracer_d;
    amrex::Vector<amrex::BCRec> m_bcrec_force;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_bcrec_force_d;

    amrex::Vector<int> m_iconserv_velocity;
    amrex::Gpu::DeviceVector<int> m_iconserv_velocity_d;
    amrex::Vector<int> m_iconserv_density;
    amrex::Gpu::DeviceVector<int> m_iconserv_density_d;
    amrex::Vector<int> m_iconserv_tracer;
    amrex::Gpu::DeviceVector<int> m_iconserv_tracer_d;

    int m_ntrac = 1;

    std::unique_ptr<DiffusionTensorOp> m_diffusion_tensor_op;
    std::unique_ptr<DiffusionScalarOp> m_diffusion_scalar_op;

    //
    // end of member variables
    //

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }
#ifdef AMREX_USE_EB
    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }
#endif

    // Number of ghost cells for field arrays.
    int nghost_state () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular()) return (m_use_godunov) ? 5 : 4;
#endif
        return (m_use_godunov) ? 3 : 2;
    }

    int nghost_force () const { return (m_use_godunov) ? 1 : 0; }

    int nghost_mac () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular()) return (m_use_godunov) ? 4 : 3;
#endif
        return (m_use_godunov) ? 1 : 0;
    }

#ifdef AMREX_USE_EB
    int nghost_eb_basic () const {
        return (m_use_godunov) ? 5 : 4;
    }

    int nghost_eb_volume () const {
        return (m_use_godunov) ? 5 : 4;
    }

    int nghost_eb_full () const {
        return (m_use_godunov) ? 3 : 2;
    }
#endif

    bool need_divtau () const {
        return !(!m_use_godunov and DiffusionType::Implicit == m_diff_type);
    }

    DiffusionTensorOp* get_diffusion_tensor_op ();
    DiffusionScalarOp* get_diffusion_scalar_op ();

    amrex::Vector<amrex::MultiFab*> get_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_vel_forces () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tra_forces () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_divtau_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_divtau_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_laps_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_laps_new () noexcept;
    //
    amrex::Vector<amrex::MultiFab const*> get_velocity_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_velocity_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_density_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_density_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tracer_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tracer_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_vel_forces_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tra_forces_const () const noexcept;

    amrex::Vector<int> const& get_velocity_iconserv () const noexcept { return m_iconserv_velocity; }
    amrex::Vector<int> const& get_density_iconserv () const noexcept { return m_iconserv_density; }
    amrex::Vector<int> const& get_tracer_iconserv () const noexcept { return m_iconserv_tracer; }
    //
    int const* get_velocity_iconserv_device_ptr () const noexcept {
        return m_iconserv_velocity_d.data(); }
    int const* get_density_iconserv_device_ptr () const noexcept {
        return m_iconserv_density_d.data(); }
    int const* get_tracer_iconserv_device_ptr () const noexcept {
        return m_iconserv_tracer_d.data(); }

    amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
    amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }
    amrex::Vector<amrex::BCRec> const& get_force_bcrec () const noexcept { return m_bcrec_force; }
    //
    amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept {
        return m_bcrec_velocity_d.data(); }
    amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept {
        return m_bcrec_density_d.data(); }
    amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept {
        return m_bcrec_tracer_d.data(); }
    amrex::BCRec const* get_force_bcrec_device_ptr () const noexcept {
        return m_bcrec_force_d.data(); }

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_projection_bc (amrex::Orientation::Side side) const noexcept;

    amrex::Vector<amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> >
    get_diffuse_tensor_bc (amrex::Orientation::Side side) const noexcept;

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_diffuse_scalar_bc (amrex::Orientation::Side side) const noexcept;

    void fillpatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillpatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillpatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);
    void fillpatch_gradp (int lev, amrex::Real time, amrex::MultiFab& gradp, int ng);
    void fillpatch_force (amrex::Real time, amrex::Vector<amrex::MultiFab*> const& force, int ng);

    void fillcoarsepatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillcoarsepatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillcoarsepatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);
    void fillcoarsepatch_gradp (int lev, amrex::Real time, amrex::MultiFab& gradp, int ng);

    void fillphysbc_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillphysbc_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillphysbc_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);

    void copy_from_new_to_old_velocity (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_velocity (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_density  (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_density  (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_tracer   (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_new_to_old_tracer   (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    //
    void copy_from_old_to_new_velocity (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_velocity (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_density  (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_density  (int lev, amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_tracer   (         amrex::IntVect const& ng = amrex::IntVect{0});
    void copy_from_old_to_new_tracer   (int lev, amrex::IntVect const& ng = amrex::IntVect{0});

    void Advance ();
    bool writeNow ();
    void AverageDown ();
    void AverageDownTo (int crse_lev);

    ///////////////////////////////////////////////////////////////////////////
    //
    // boundary_conditions
    //
    ///////////////////////////////////////////////////////////////////////////

    void init_bcs ();

    ///////////////////////////////////////////////////////////////////////////
    //
    // convection
    //
    ///////////////////////////////////////////////////////////////////////////

    void init_advection ();

    void apply_MAC_projection (amrex::Vector<amrex::MultiFab*> const& u_mac,
                               amrex::Vector<amrex::MultiFab*> const& v_mac,
                               amrex::Vector<amrex::MultiFab*> const& w_mac,
                               amrex::Vector<amrex::MultiFab const*> const& density,
                               amrex::Real time);

    ///////////////////////////////////////////////////////////////////////////
    //
    // embedded_boundaries
    //
    ///////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_EB
    void MakeEBGeometry ();
    void make_eb_annulus ();
    void make_eb_box ();
    void make_eb_cylinder ();
    void make_eb_twocylinders ();
    void make_eb_regular ();
    void make_eb_sphere ();
    void make_eb_spherecube ();
    void make_eb_tuscan ();
#endif

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
