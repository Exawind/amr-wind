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

public: // for cuda

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
                                amrex::MultiFab const& vel);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Godunov Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

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

    void predict_godunov_on_box (int lev, amrex::Box const& bx, int ncomp,
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
                   amrex::FabFactory<amrex::FArrayBox> const& fact,
                   int ntrac, int ng_state, int ng_force);
        // cell-centered multifabs
        amrex::MultiFab velocity;
        amrex::MultiFab velocity_o;
        amrex::MultiFab density;
        amrex::MultiFab density_o;
        amrex::MultiFab tracer;
        amrex::MultiFab tracer_o;
        amrex::MultiFab gp; // pressure gradient
        amrex::MultiFab vel_forces;
        amrex::MultiFab tra_forces;
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

    amrex::FabFactory<amrex::FArrayBox> const&
    Factory (int lev) const noexcept { return *m_factory[lev]; }
#ifdef AMREX_USE_EB
    amrex::EBFArrayBoxFactory const&
    EBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
    }
#endif

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

    int ntrac = 1;

    // Number of ghost cells for field arrays.
    int nghost_state () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular()) return (use_godunov) ? 5 : 4;
#endif
        return (use_godunov) ? 3 : 2;
    }

    int nghost_force () const { return (use_godunov) ? 1 : 0; }

    int nghost_mac () const {
#ifdef AMREX_USE_EB
        if (!EBFactory(0).isAllRegular()) return 3;
#endif
        return 0;
    }

#ifdef AMREX_USE_EB
    int nghost_eb_basic () const {
        return (use_godunov) ? 5 : 4;
    }

    int nghost_eb_volume () const {
        return (use_godunov) ? 5 : 4;
    }

    int nghost_eb_full () const {
        return (use_godunov) ? 3 : 2;
    }
#endif

    enum struct DiffusionType {
        Invalid, Explicit, Crank_Nicolson, Implicit
    };
    DiffusionType m_diff_type = DiffusionType::Implicit;
    std::unique_ptr<DiffusionTensorOp> diffusion_tensor_op;
    std::unique_ptr<DiffusionScalarOp> diffusion_scalar_op;

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
    //
    amrex::Vector<amrex::MultiFab const*> get_velocity_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_velocity_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_density_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_density_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tracer_old_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tracer_new_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_vel_forces_const () const noexcept;
    amrex::Vector<amrex::MultiFab const*> get_tra_forces_const () const noexcept;

    amrex::Vector<amrex::BCRec> const& get_velocity_bcrec () const noexcept { return m_bcrec_velocity; }
    amrex::Vector<amrex::BCRec> const& get_density_bcrec () const noexcept { return m_bcrec_density; }
    amrex::Vector<amrex::BCRec> const& get_tracer_bcrec () const noexcept { return m_bcrec_tracer; }

    amrex::BCRec const* get_velocity_bcrec_device_ptr () const noexcept {
        return m_bcrec_velocity_d.data(); }
    amrex::BCRec const* get_density_bcrec_device_ptr () const noexcept {
        return m_bcrec_density_d.data(); }
    amrex::BCRec const* get_tracer_bcrec_device_ptr () const noexcept {
        return m_bcrec_tracer_d.data(); }

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_projection_bc (amrex::Orientation::Side side) const noexcept;

    amrex::Vector<amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> >
    get_diffuse_tensor_bc (amrex::Orientation::Side side) const noexcept;

    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
    get_diffuse_scalar_bc (amrex::Orientation::Side side) const noexcept;

    void fillpatch_velocity (int lev, amrex::Real time, amrex::MultiFab& vel, int ng);
    void fillpatch_density (int lev, amrex::Real time, amrex::MultiFab& density, int ng);
    void fillpatch_tracer (int lev, amrex::Real time, amrex::MultiFab& tracer, int ng);

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
                                amrex::Box const& domain,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

    void init_double_shear_layer (amrex::Box const& vbx, amrex::Box const& gbx,
                                  amrex::Array4<amrex::Real> const& p,
                                  amrex::Array4<amrex::Real> const& vel,
                                  amrex::Array4<amrex::Real> const& denisty,
                                  amrex::Array4<amrex::Real> const& tracer,
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

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Initialization
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Read inputs file using ParmParse
    void ReadParameters ();
    void ReadIOParameters();
    void ReadRheologyParameters();

    // Resize arrays to fit (up to) max_level + 1 AMR levels
    void ResizeArrays();

    // Initialization
    void InitialProjection();
    void InitialIterations();

    // Member variables for initial conditions
    int probtype = 0;
    amrex::Real ic_u = 0.0;
    amrex::Real ic_v = 0.0;
    amrex::Real ic_w = 0.0;
    amrex::Real ic_p = 0.0;

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

    void compute_forces (amrex::Vector<amrex::MultiFab*> const& vel_forces,
                         amrex::Vector<amrex::MultiFab*> const& tra_forces,
                         amrex::Vector<amrex::MultiFab const*> const& density,
                         amrex::Vector<amrex::MultiFab const*> const& tracer);

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
                                  amrex::Vector<amrex::MultiFab const*> const& vel,
                                  amrex::Vector<amrex::MultiFab const*> const& density,
                                  amrex::Vector<amrex::MultiFab const*> const& tracer,
                                  amrex::Vector<amrex::MultiFab const*> const& vel_forces,
                                  amrex::Vector<amrex::MultiFab const*> const& tra_forces,
                                  amrex::Real time);

    void apply_MAC_projection (amrex::Vector<amrex::MultiFab>& u_mac,
                               amrex::Vector<amrex::MultiFab>& v_mac,
                               amrex::Vector<amrex::MultiFab>& w_mac,
                               amrex::Vector<amrex::MultiFab const*> const& density,
                               amrex::Real time);

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

    void make_eb_annulus();
    void make_eb_box();
    void make_eb_cylinder();
    void make_eb_twocylinders();
    void make_eb_regular();
    void make_eb_sphere();
    void make_eb_spherecube();
    void make_eb_tuscan();
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
    amrex::Real init_shrink = 0.1;

    // Initial projection / iterations
    bool do_initial_proj    = true;
    int  initial_iterations = 3;

    // Use Boussinesq approximation for buoyancy?
    bool use_boussinesq = false;

    // AMR / refinement settings
    int refine_cutcells = 1;
    int regrid_int = -1;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Physics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

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

    bool use_godunov         = false;
    bool use_forces_in_trans = false;
    int  diffusion_type      = 2;

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

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diagnostics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void PrintMaxValues(amrex::Real time);
    void PrintMaxVel(int lev);
    void PrintMaxGp(int lev);
    void CheckForNans(int lev);

    bool writeNow();

    void AverageDown();
    void AverageDownTo(int crse_lev);
};

#endif
