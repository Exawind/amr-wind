#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_NodalProjector.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <eb_if.H>
#endif

#include <DiffusionOp.H>
#include <incflo_bclist.H>

class incflo : public amrex::AmrCore
{
public:

    incflo ();
    virtual ~incflo();

    // Initialize multilevel AMR data
    void InitData();

#ifdef AMREX_USE_EB
    void WriteMyEBSurface();
#endif

    // Evolve solution to final time through repeated calls to Advance()
    void Evolve();

    // Tag cells for refinement
    virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping
    // Only used during initialization
    virtual void MakeNewLevelFromScratch(int lev, amrex::Real time,
                                         const amrex::BoxArray& new_grids,
                                         const amrex::DistributionMapping& new_dmap) override;

    // Make a new level using provided BoxArray and DistributionMapping,
    // and fill with interpolated coarse level data
    virtual void MakeNewLevelFromCoarse(int lev, amrex::Real time,
                                        const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping,
    // and fill with existing fine and coarse data
    virtual void RemakeLevel(int lev, amrex::Real time,
                             const amrex::BoxArray& ba,
                             const amrex::DistributionMapping& dm) override;

    // Delete level data
    virtual void ClearLevel(int lev) override;

    // Set up nodal and diffusion solvers
    void incflo_init_solvers();
    void incflo_setup_solvers();

    const int* get_bc_ilo_ptr(int lev){ return bc_ilo[lev]->dataPtr(); }
    const int* get_bc_ihi_ptr(int lev){ return bc_ihi[lev]->dataPtr(); }
    const int* get_bc_jlo_ptr(int lev){ return bc_jlo[lev]->dataPtr(); }
    const int* get_bc_jhi_ptr(int lev){ return bc_jhi[lev]->dataPtr(); }
    const int* get_bc_klo_ptr(int lev){ return bc_klo[lev]->dataPtr(); }
    const int* get_bc_khi_ptr(int lev){ return bc_khi[lev]->dataPtr(); }

    int get_nghost(){ return nghost; }
    int get_probtype(){ return probtype; }
    void GetInputBCs();

    void incflo_set_velocity_bcs (amrex::Real time,
                                  amrex::Vector< std::unique_ptr<amrex::MultiFab> >& vel);

    void incflo_set_density_bcs (amrex::Real time,
                                 amrex::Vector< std::unique_ptr<amrex::MultiFab> >& density);

    void incflo_set_tracer_bcs (amrex::Real time,
                                amrex::Vector< std::unique_ptr<amrex::MultiFab> >& trac);

    void set_density_bcs(amrex::Real time,
                         const int lev,
                         amrex::FArrayBox& scal_fab,
                         const amrex::Box& domain);

    void set_tracer_bcs(amrex::Real time,
                        const int lev,
                        amrex::FArrayBox& scal_fab,
                        const int comp,
                        const int ncomp,
                        const amrex::Box& domain);

    void set_velocity_bcs(amrex::Real time,
                          const int lev,
                          amrex::FArrayBox& vel_fab,
                          const amrex::Box& domain) const;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeStrainrate(amrex::Real time_in);

    void ComputeViscosity( amrex::Vector<std::unique_ptr<amrex::MultiFab>>& eta_out,
                           const amrex::Real time_in);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Derived quantities and diagnostics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeVorticity(amrex::Real time_in);
    double ComputeKineticEnergy();

    amrex::Real volWgtSum (int lev, const amrex::MultiFab & mf, int comp, bool local=false);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // MOL Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void predict_vels_on_faces (int lev, amrex::Real time, amrex::MultiFab& u_mac,
                                amrex::MultiFab& v_mac, amrex::MultiFab& w_mac,
                                amrex::MultiFab& vel);

    void incflo_predict_vels_on_faces(int lev, amrex::Real time,
                                      amrex::Vector< std::unique_ptr<amrex::MultiFab> >& vel);

    void incflo_compute_slopes(int lev, amrex::Real time, amrex::MultiFab& Sborder,
                               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& xslopes_in,
                               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& yslopes_in,
                               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& zslopes_in,
                               int slopes_comp, int ncomp);

    void set_MAC_velocity_bcs (int lev,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u_mac,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v_mac,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w_mac,
                               amrex::Real time);

    void incflo_compute_fluxes(const int lev,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& a_fx,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& a_fy,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& a_fz,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& state,
                               const int state_comp, const int n_comp,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& xslopes,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& yslopes,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& zslopes,
                               const int slopes_comp,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u_mac,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v_mac,
                               amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w_mac);

    void incflo_compute_fluxes_on_box(const int lev, amrex::Box& bx,
                                      amrex::FArrayBox& a_fx,
                                      amrex::FArrayBox& a_fy,
                                      amrex::FArrayBox& a_fz,
                                      const amrex::FArrayBox& state,
                                      const int state_comp, const int n_comp,
                                      const amrex::FArrayBox& xslopes,
                                      const amrex::FArrayBox& yslopes,
                                      const amrex::FArrayBox& zslopes,
                                      const int slopes_comp,
                                      const amrex::FArrayBox& u_mac,
                                      const amrex::FArrayBox& v_mac,
                                      const amrex::FArrayBox& w_mac);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Godunov Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void incflo_predict_godunov(int lev, amrex::Real time,
                                amrex::Vector< std::unique_ptr<amrex::MultiFab> >& vel);

    void incflo_predict_godunov_on_box (const int lev, amrex::Box& bx,
                                        const amrex::Array4<      amrex::Real> &a_vel,
                                        const amrex::Array4<const amrex::Real> &u_face,
                                        const amrex::Array4<const amrex::Real> &v_face,
                                        const amrex::Array4<const amrex::Real> &w_face,
                                        const amrex::Array4<      amrex::Real> &tf,
                                        const amrex::Array4<      amrex::Real> &divu_cc,
                                        const amrex::Gpu::ManagedVector<amrex::BCRec> &BCs);

    void incflo_godunov_fluxes_on_box(const int lev, amrex::Box& bx,
                                      const amrex::Array4<amrex::Real> &a_fx,
                                      const amrex::Array4<amrex::Real> &a_fy,
                                      const amrex::Array4<amrex::Real> &a_fz,
                                      const amrex::Array4<amrex::Real> &tforces,
                                      const amrex::Array4<amrex::Real> &divu,
                                      const amrex::Array4<amrex::Real> &s,
                                      const int state_comp, const int n_comp,
                                      const amrex::Array4<const amrex::Real> &umac,
                                      const amrex::Array4<const amrex::Real> &vmac,
                                      const amrex::Array4<const amrex::Real> &wmac,
                                      const amrex::Gpu::ManagedVector<amrex::BCRec> &BCs);

#ifdef AMREX_USE_EB
    void incflo_compute_eb_fluxes_on_box(const int lev, amrex::Box& bx,
                                         amrex::FArrayBox& a_fx,
                                         amrex::FArrayBox& a_fy,
                                         amrex::FArrayBox& a_fz,
                                         const amrex::FArrayBox& state,
                                         const int state_comp, const int n_comp,
                                         const amrex::FArrayBox& xslopes,
                                         const amrex::FArrayBox& yslopes,
                                         const amrex::FArrayBox& zslopes,
                                         const int slopes_comp,
                                         const amrex::FArrayBox& ep_u_mac,
                                         const amrex::FArrayBox& ep_v_mac,
                                         const amrex::FArrayBox& ep_w_mac,
                                         const amrex::FArrayBox& ax,
                                         const amrex::FArrayBox& ay,
                                         const amrex::FArrayBox& az,
                                         const amrex::FArrayBox& fc_x,
                                         const amrex::FArrayBox& fc_y,
                                         const amrex::FArrayBox& fc_z,
                                         const amrex::FArrayBox& volfrac,
                                         const amrex::FArrayBox& bndrycent,
                                         const amrex::IArrayBox& cc_mask,
                                         const amrex::EBCellFlagFab& flags);
#endif

private:

    // Be verbose?
    int incflo_verbose = 0;

    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> t_new;

    // Times
    amrex::Real cur_time     = 0.0;
    amrex::Real dt           = -1.0;
    amrex::Real prev_dt      = -1.0;
    amrex::Real prev_prev_dt = -1.0;

    struct LevelData {
        LevelData () = default;
        LevelData (amrex::BoxArray const& ba,
                   amrex::DistributionMapping const& dm,
                   amrex::FabFactory<amrex::FArrayBox> const& fact);
        // cell-centered multifabs
        amrex::MultiFab velocity;
        amrex::MultiFab velocity_o;
        amrex::MultiFab density;
        amrex::MultiFab density_o;
        amrex::MultiFab tracer;
        amrex::MultiFab tracer_o;
        amrex::MultiFab gp; // pressure gradient
        amrex::MultiFab eta; // xxxx seems that we do not need it for "newtonian" fluid model.
        // nodal multifabs
        amrex::MultiFab p;
        //
        amrex::MultiFab conv_velocity;
        amrex::MultiFab conv_velocity_o;
        amrex::MultiFab conv_density;
        amrex::MultiFab conv_density_o;
        amrex::MultiFab conv_tracer;
        amrex::MultiFab conv_tracer_o;
    };

    amrex::Vector<std::unique_ptr<LevelData> > m_leveldata;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

    enum struct BC {
        pressure_inflow, pressure_outflow, mass_inflow, no_slip_wall, periodic, undefined
    };

    amrex::GpuArray<BC                         , AMREX_SPACEDIM*2> m_bc_type;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_pressure;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_density;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
                                               , AMREX_SPACEDIM*2> m_bc_velocity;
    amrex::GpuArray<amrex::Vector<amrex::Real> , AMREX_SPACEDIM*2> m_bc_tracer;
    amrex::GpuArray<amrex::Real                , AMREX_SPACEDIM*2> m_bc_location;

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

    static int ntrac;

    // Number of ghost nodes for field arrays.
#ifdef AMREX_USE_EB
    static constexpr int nghost            = 4;
    static constexpr int nghost_for_slopes = 4;
    static constexpr int nghost_for_bcs    = 4;
#else
    static constexpr int nghost            = 2;
    static constexpr int nghost_for_slopes = 2;
    static constexpr int nghost_for_bcs    = 2;
#endif

#ifdef AMREX_USE_EB
    static constexpr int m_eb_basic_grow_cells = nghost;
    static constexpr int m_eb_volume_grow_cells = nghost;
    static constexpr int m_eb_full_grow_cells = 2;
#endif

    std::unique_ptr<amrex::NodalProjector> nodal_projector;

    enum struct DiffusionType {
        Invalid, Explicit, Crank_Nicolson, Implicit
    };
    DiffusionType m_diff_type = DiffusionType::Implicit;
    std::unique_ptr<DiffusionOp> diffusion_op;

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }
#ifdef AMREX_USE_EB
    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }
#endif

    amrex::Vector<amrex::MultiFab*> get_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_tracer_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_velocity_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_density_new () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_old () noexcept;
    amrex::Vector<amrex::MultiFab*> get_conv_tracer_new () noexcept;

    amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
    amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }

    amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept {
        return m_bcrec_velocity_d.data(); }
    amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept {
        return m_bcrec_density_d.data(); }
    amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept {
        return m_bcrec_tracer_d.data(); }

    std::array<amrex::LinOpBCType,AMREX_SPACEDIM> get_projection_bclo () const noexcept;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM> get_projection_bchi () const noexcept;

    void fillpatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillpatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillpatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);

public: // for cuda

    void init_bcs ();
    void set_background_pressure ();
    void prob_init_fluid (int lev);
    void prob_set_inflow_velocity (int grid_id, amrex::Orientation ori, amrex::Box const& bx,
                                   amrex::Array4<amrex::Real> const& v, int lev, amrex::Real time);

    void init_plane_poiseuille (amrex::Box const& vbx, amrex::Box const& gbx,
                                amrex::Array4<amrex::Real> const& p,
                                amrex::Array4<amrex::Real> const& vel,
                                amrex::Array4<amrex::Real> const& denisty,
                                amrex::Array4<amrex::Real> const& tracer,
                                amrex::Array4<amrex::Real> const& eta,
                                amrex::Box const& domain,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void set_inflow_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int nghost);

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
                                   amrex::Array4<amrex::Real const> const& fcz);
#endif

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_u;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_v;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_w;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_r;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_t;
    mutable amrex::Gpu::ManagedVector<amrex::Real> m_bc_p;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Initialization
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Set boundary conditions from inputs file
    void SetInputBCs(const std::string bcID, const int index, const int cyclic, const amrex::Real domloc);

    // Read inputs file using ParmParse
    void ReadParameters();
    void ReadIOParameters();
    void ReadRheologyParameters();

    // Resize arrays to fit (up to) max_level + 1 AMR levels
    void ResizeArrays();

    // Post-initialization: set BCs, apply ICs, initial velocity projection, pressure iterations
    void PostInit(int restart_flag);
    void SetBCTypes();
    void InitFluid();
    void SetBackgroundPressure();
    void InitialProjection();
    void InitialIterations();
    void setup_level_mask();

    // Member variables for initial conditions
    int probtype = 0;
    amrex::Real ic_u = 0.0;  // xxxxx Do we really need it?
    amrex::Real ic_v = 0.0;
    amrex::Real ic_w = 0.0;
    amrex::Real ic_p = 0.0;

    // LinOpBCType Definitions from amrex/Src/Boundary/AMReX_LO_BCTYPES.H -- these for projections
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  ppe_lobc;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  ppe_hibc;

    // LinOpBCType Definitions from amrex/Src/Boundary/AMReX_LO_BCTYPES.H -- these for velocity diffusion solves
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  vel_diff_lobc;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  vel_diff_hibc;

    // LinOpBCType Definitions from amrex/Src/Boundary/AMReX_LO_BCTYPES.H -- these for scalar diffusion solves
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  scal_diff_lobc;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  scal_diff_hibc;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Simulate!
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void Advance();
    void ComputeDt(int initialisation, bool explicit_diffusion);
    bool SteadyStateReached();
    void ApplyPredictor(bool incremental_projection = false);
    void ApplyCorrector();
    void ApplyProjection(amrex::Real time, amrex::Real scaling_factor, bool incremental);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

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
                                       amrex::Array4<amrex::Real> const& qface);

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

    void compute_convective_term (amrex::Box const& bx, int lev, amrex::MFIter const& mfi,
                                  amrex::Array4<amrex::Real> const& dvdt,
                                  amrex::Array4<amrex::Real> const& drdt,
                                  amrex::Array4<amrex::Real> const& dtdt,
                                  amrex::Array4<amrex::Real const> const& vel,
                                  amrex::Array4<amrex::Real const> const& rho,
                                  amrex::Array4<amrex::Real const> const& tra,
                                  amrex::Array4<amrex::Real const> const& umac,
                                  amrex::Array4<amrex::Real const> const& vmac,
                                  amrex::Array4<amrex::Real const> const& wmac);

    void compute_convective_term (amrex::Vector<amrex::MultiFab*> const& conv_u,
                                  amrex::Vector<amrex::MultiFab*> const& conv_r,
                                  amrex::Vector<amrex::MultiFab*> const& conv_t,
                                  amrex::Vector<amrex::MultiFab*> const& vel,
                                  amrex::Vector<amrex::MultiFab*> const& density,
                                  amrex::Vector<amrex::MultiFab*> const& tracer,
                                  amrex::Real time);

    void apply_MAC_projection (amrex::Vector<amrex::MultiFab>& u_mac,
                               amrex::Vector<amrex::MultiFab>& v_mac,
                               amrex::Vector<amrex::MultiFab>& w_mac,
                               amrex::Vector<amrex::MultiFab*> const& density,
                               amrex::Real time);

     void incflo_compute_convective_term(amrex::Vector< std::unique_ptr<amrex::MultiFab> >& conv_u,
                                         amrex::Vector< std::unique_ptr<amrex::MultiFab> >& conv_r,
                                         amrex::Vector< std::unique_ptr<amrex::MultiFab> >& conv_t,
                                         amrex::Vector< std::unique_ptr<amrex::MultiFab> >& vel,
                                         amrex::Vector< std::unique_ptr<amrex::MultiFab> >& density_in,
                                         amrex::Vector< std::unique_ptr<amrex::MultiFab> >& tracer_in,
                                         amrex::Real time);

    void apply_MAC_projection ( amrex::Vector< std::unique_ptr<amrex::MultiFab> >& u_mac,
                                amrex::Vector< std::unique_ptr<amrex::MultiFab> >& v_mac,
                                amrex::Vector< std::unique_ptr<amrex::MultiFab> >& w_mac,
                                amrex::Vector< std::unique_ptr<amrex::MultiFab> >& density,
                                amrex::Real time);

     void incflo_divergence_plus_redist(const int lev, amrex::Vector< std::unique_ptr<amrex::MultiFab> >& conv,
                                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& fx,
                                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& fy,
                                        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& fz,
                                        int num_comp);

     // Verbosity
     int  mac_mg_verbose = 0;
     int  mac_mg_cg_verbose = 0;

     // Control MLMG behavior
     int  mac_mg_maxiter = 200;
     int  mac_mg_cg_maxiter = 200;
     amrex::Real mac_mg_rtol = 1.0e-11;
     amrex::Real mac_mg_atol = 1.0e-14;
     std::string mac_bottom_solver_type;

     // Max coarsening level
     int  mac_mg_max_coarsening_level = 100;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Nodal Projection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeDivU(amrex::Real time_in);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Derived quantities
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeDrag();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Embedded Boundaries
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_EB
    void MakeEBGeometry();
    bool UpdateEBFactory(int a_lev);

    void make_eb_annulus();
    void make_eb_box();
    void make_eb_cylinder();
    void make_eb_twocylinders();
    void make_eb_regular();
    void make_eb_sphere();
    void make_eb_spherecube();
    void make_eb_tuscan();

    const amrex::EB2::Level* eb_level;
#endif

#if 0
    // If AMREX_USE_EB = TRUE this will be be an EBFArrayBoxFactory,
    //    but otherwise it will just be a FabFactory that doesn't know anything about EB
    //    amrex::Vector<std::unique_ptr<FabFactory<FArrayBox>>> ebfactory;
#else
    //  amrex::Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;
#endif

#ifdef AMREX_USE_EB
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>> ebfactory;
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Runtime parameters
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Time step counter
    int nstep         = -1;

    // Stop simulation if cur_time reaches stop_time OR nstep reaches max_step
    // OR steady_state = true AND steady_state_tol is reached
    amrex::Real stop_time = -1.0;
    int max_step = -1;
    bool steady_state = false;
    amrex::Real steady_state_tol = 1.0e-5;

    // Options to control time stepping
    amrex::Real cfl = 0.5;
    amrex::Real fixed_dt = -1.;

    // Initial projection / iterations
    bool do_initial_proj    = true;
    int  initial_iterations = 3;

    // Use Boussinesq approximation for buoyancy?
    bool use_boussinesq = false;

    // AMR / refinement settings
    int refine_cutcells = 1;
    int regrid_int = -1;

    BCList bc_list;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Physics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    amrex::Periodicity p0_periodicity;
    amrex::Real p000 = 0.0; // background pressure at node (0,0,0).
    amrex::Vector<amrex::Real> gp0{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> gravity{{0.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> delp{{0.0, 0.0, 0.0}};

    // Scalar diffusive coefficient
    amrex::Vector<amrex::Real> mu_s;

    // Density (if constant)
    amrex::Real ro_0 = 1.0;

    bool constant_density = true;
    bool advect_tracer    = false;
    bool test_tracer_conservation = false;
    bool use_godunov      = false;

    // Fluid properties
    std::string fluid_model;
    amrex::Real mu = 1.0;
    amrex::Real n_0 = 0.0;
    amrex::Real tau_0 = 0.0;
    amrex::Real papa_reg = 0.0;
    amrex::Real eta_0 = 0.0;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Input / Output
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void WriteHeader(const std::string& name, bool is_checkpoint) const;
    void WriteJobInfo(const std::string& dir) const;
    void WriteCheckPointFile() const;
    void WritePlotFile();
    void ReadCheckpointFile();

    // Member variables for I/O
    int plot_int = -1;

    // Dump plotfiles at as close as possible to the designated period *without* changing dt
    amrex::Real plot_per_approx = -1.0;
    
    // Dump plotfiles at exactcly the designated period by changing dt
    amrex::Real plot_per_exact  = -1.0;

    int last_plt = -1;
    std::string plot_file{"plt"};

    int check_int = -1;
    int last_chk = -1;
    int KE_int = -1;
    std::string check_file{"chk"};
    std::string restart_file{""};

    // Flags for saving fluid data in plot files
    int plt_velx        = 1;
    int plt_vely        = 1;
    int plt_velz        = 1;
    int plt_gpx         = 1;
    int plt_gpy         = 1;
    int plt_gpz         = 1;
    int plt_rho         = 1;
    int plt_tracer      = 1;
    int plt_p           = 0;
    int plt_eta         = 0;
    int plt_vort        = 1;
    int plt_strainrate  = 0;
    int plt_stress      = 0;
    int plt_divu        = 0;
    int plt_vfrac       = 1;

    // Total number of variables to write in plot file
    int pltVarCount = 0;

    // amrex::Vector variables to simplify IO -- the vectors are the same for plt/chk
    amrex::Vector<const amrex::Vector<std::unique_ptr<amrex::MultiFab>>*> vectorVars;
    amrex::Vector<std::string> vecVarsName = {"velx", "vely", "velz", "gpx", "gpy", "gpz"};

    // Scalar variables at cell centers that need to be written to checkfile.
    amrex::Vector<const amrex::Vector<std::unique_ptr<amrex::MultiFab>>*> chkscalarVars = {&p, &density, &tracer, &eta};
    amrex::Vector<std::string> chkscaVarsName = {"p", "density", "tracer", "eta"};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diagnostics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    amrex::Real Norm(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mf, int lev, int comp, int norm_type);
    void PrintMaxValues(amrex::Real time);
    void PrintMaxVel(int lev);
    void PrintMaxGp(int lev);
    void CheckForNans(int lev);

    const amrex::Real covered_val = 1.0e40;
    const amrex::Real boundary_val = 1.0e80;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Data containers
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Boundary conditions
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_ilo;
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_ihi;
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_jlo;
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_jhi;
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_klo;
    amrex::Vector<std::unique_ptr<amrex::IArrayBox>> bc_khi;

    // Primary variables
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> density;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> density_o;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> tracer;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> tracer_o;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_o;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> p;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> p0;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> gp;

    // xxxxx could we make some amrex::MultiFabs here local to functions

    // Derived variables TODO: to save memory, would only need to have 2 temporary variables
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> eta;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> eta_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> strainrate;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vort;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> drag;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> divu;

    // Helper variables
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_u_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_r_old;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_t;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> conv_t_old;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> divtau;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> divtau_old;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> laps;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> laps_old;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> xslopes_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> yslopes_u;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> zslopes_u;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> xslopes_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> yslopes_r;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> zslopes_r;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> xslopes_t;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> yslopes_t;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> zslopes_t;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_u_mac;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_v_mac;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_w_mac;

    // integer amrex::MultiFab specifying if a finer level exists above
    // value is 1 if a finer level is not above and 0 otherwise
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> level_mask;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Not yet classified
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    bool writeNow();

    void FillPatchVel    (int lev, amrex::Real time, amrex::MultiFab& mf);
    void FillPatchDensity(int lev, amrex::Real time, amrex::MultiFab& mf);
    void FillPatchScalar (int lev, amrex::Real time, amrex::MultiFab& mf);

    void GetDataVel    (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);
    void GetDataDensity(int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);
    void GetDataScalar (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);

    void AverageDown();
    void AverageDownTo(int crse_lev);

    void AllocateArrays(int lev);
    void RegridArrays(int lev);
    void MakeBCArrays();

};

#endif
