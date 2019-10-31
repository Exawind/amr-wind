#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_NodalProjector.H>

#include <eb_if.H>
#include <DiffusionOp.H>
#include <incflo_bclist.H>


enum DiffusionType
{
    Invalid=-1,
    Explicit,
    Crank_Nicolson,
    Implicit
};

class incflo : public AmrCore
{
public:
    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Public member functions are defined directly in incflo.cpp.
    // They include the functions called directly in main.cpp,
    // and the 5 pure virtual functions in the AmrCore class, which must be implemented.
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Functions called directly in main()
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Constructor
    //  - reads in parameters from inputs file
    //  - sizes multilevel arrays and datastructures
    //  - makes EB geometry
	incflo();
	~incflo();

    // Initialize multilevel AMR data
    void InitData();
    BoxArray MakeBaseGrids () const;
    void ChopGrids (const Box& domain, BoxArray& ba, int target_size) const;

    void WriteMyEBSurface();

    // Evolve solution to final time through repeated calls to Advance()
    void Evolve();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Pure virtual functions from AmrCore class
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Tag cells for refinement
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping
    // Only used during initialization
    void MakeNewLevelFromScratch(int lev, Real time,
                                 const BoxArray& new_grids,
                                 const DistributionMapping& new_dmap) override;

    // Make a new level using provided BoxArray and DistributionMapping,
    // and fill with interpolated coarse level data
    void MakeNewLevelFromCoarse(int lev, Real time,
                                const BoxArray& ba,
                                const DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping,
    // and fill with existing fine and coarse data
    void RemakeLevel(int lev, Real time,
                             const BoxArray& ba,
                             const DistributionMapping& dm) override;

    // Delete level data
    void ClearLevel(int lev) override;

    // Set up nodal and diffusion solvers
    void incflo_init_solvers();
    void incflo_setup_solvers();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Some getters (TODO: find better way to do fillpatching)
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    const Geometry get_geom_ref(int lev){ return geom[lev]; }

    const int* get_bc_ilo_ptr(int lev){ return bc_ilo[lev]->dataPtr(); }
    const int* get_bc_ihi_ptr(int lev){ return bc_ihi[lev]->dataPtr(); }
    const int* get_bc_jlo_ptr(int lev){ return bc_jlo[lev]->dataPtr(); }
    const int* get_bc_jhi_ptr(int lev){ return bc_jhi[lev]->dataPtr(); }
    const int* get_bc_klo_ptr(int lev){ return bc_klo[lev]->dataPtr(); }
    const int* get_bc_khi_ptr(int lev){ return bc_khi[lev]->dataPtr(); }

    int get_nghost(){ return nghost; }
    int get_probtype(){ return probtype; }
    void GetInputBCs();

    void incflo_set_velocity_bcs (Real time,
                                  Vector< std::unique_ptr<MultiFab> >& vel,
                                  int extrap_dir_bcs) const;

    void incflo_set_density_bcs (Real time,
                                 Vector< std::unique_ptr<MultiFab> >& density);

    void incflo_set_tracer_bcs (Real time,
                                Vector< std::unique_ptr<MultiFab> >& trac);

    void set_density_bcs(Real time,
                         const int lev,
                         FArrayBox& scal_fab,
                         const Box& domain);

    void set_tracer_bcs(Real time,
                        const int lev,
                        FArrayBox& scal_fab,
                        const int comp,
                        const int ncomp,
                        const Box& domain);

    void set_velocity_bcs(Real time,
                          const int lev,
                          FArrayBox& vel_fab,
                          const Box& domain,
                          const int* extrap_dir_bcs) const;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeStrainrate(Real time_in);

    void ComputeViscosity( Vector<std::unique_ptr<MultiFab>>& eta_out,
                           const Real time_in);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Derived quantities
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeVorticity(Real time_in);
    double ComputeKineticEnergy();


    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void incflo_predict_vels_on_faces(int lev, Real time,
                                      Vector< std::unique_ptr<MultiFab> >& vel);

    void incflo_compute_slopes(int lev, Real time, MultiFab& Sborder,
                               Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                               int slopes_comp, int ncomp);

    void set_MAC_velocity_bcs (int lev,
                               Vector< std::unique_ptr<MultiFab> >& mac_rhs,
                               Vector< std::unique_ptr<MultiFab> >& u_mac,
                               Vector< std::unique_ptr<MultiFab> >& v_mac,
                               Vector< std::unique_ptr<MultiFab> >& w_mac,
                               Real time);

    void incflo_compute_fluxes(const int lev,
                               Vector< std::unique_ptr<MultiFab> >& a_fx,
                               Vector< std::unique_ptr<MultiFab> >& a_fy,
                               Vector< std::unique_ptr<MultiFab> >& a_fz,
                               Vector< std::unique_ptr<MultiFab> >& state,
                               const int state_comp, const int n_comp,
                               Vector< std::unique_ptr<MultiFab> >& xslopes,
                               Vector< std::unique_ptr<MultiFab> >& yslopes,
                               Vector< std::unique_ptr<MultiFab> >& zslopes,
                               const int slopes_comp,
                               Vector< std::unique_ptr<MultiFab> >& u_mac,
                               Vector< std::unique_ptr<MultiFab> >& v_mac,
                               Vector< std::unique_ptr<MultiFab> >& w_mac);

    void incflo_compute_ugradu(const int lev, Box& bx,
                               Vector< std::unique_ptr<MultiFab> >& a_fx,
                               Vector< std::unique_ptr<MultiFab> >& a_fy,
                               Vector< std::unique_ptr<MultiFab> >& a_fz,
                               Vector< std::unique_ptr<MultiFab> >& state,
                               const int state_comp, const int n_comp,
                               Vector< std::unique_ptr<MultiFab> >& xslopes,
                               Vector< std::unique_ptr<MultiFab> >& yslopes,
                               Vector< std::unique_ptr<MultiFab> >& zslopes,
                               const int slopes_comp,
                               Vector< std::unique_ptr<MultiFab> >& u_mac,
                               Vector< std::unique_ptr<MultiFab> >& v_mac,
                               Vector< std::unique_ptr<MultiFab> >& w_mac,
                               MFIter* mfi, Box& domain);

    void incflo_compute_ugradu_eb(const int lev, Box& bx,
                                  Vector< std::unique_ptr<MultiFab> >& a_fx,
                                  Vector< std::unique_ptr<MultiFab> >& a_fy,
                                  Vector< std::unique_ptr<MultiFab> >& a_fz,
                                  Vector< std::unique_ptr<MultiFab> >& state,
                                  const int state_comp, const int n_comp,
                                  Vector< std::unique_ptr<MultiFab> >& xslopes,
                                  Vector< std::unique_ptr<MultiFab> >& yslopes,
                                  Vector< std::unique_ptr<MultiFab> >& zslopes,
                                  const int slopes_comp,
                                  Vector< std::unique_ptr<MultiFab> >& u_mac,
                                  Vector< std::unique_ptr<MultiFab> >& v_mac,
                                  Vector< std::unique_ptr<MultiFab> >& w_mac,
                                  MFIter* mfi,
                                  Array< const MultiCutFab*,AMREX_SPACEDIM>& areafrac,
                                  Array< const MultiCutFab*,AMREX_SPACEDIM>& facecent,
                                  const amrex::MultiFab* volfrac,
                                  const amrex::MultiCutFab* bndrycent,
                                  const amrex::iMultiFab* cc_mask,
                                  Box& domain, const EBCellFlagFab& flags);

private:

     static int nlev;
     static int ntrac;

     // Unit vectors in Cartesian space
     static IntVect e_x;
     static IntVect e_y;
     static IntVect e_z;

     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_u;
     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_v;
     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_w;
     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_r;
     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_t;
     mutable amrex::Cuda::ManagedVector<amrex::Real> m_bc_p;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Initialization
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Set boundary conditions from inputs file
    void SetInputBCs(const std::string bcID, const int index, const int cyclic, const Real domloc);

    // Read inputs file using ParmParse
    void ReadParameters();

    // Resize arrays to fit (up to) max_level + 1 AMR levels
	void ResizeArrays();

    // Post-initialization: set BCs, apply ICs, initial velocity projection, pressure iterations, level_mask
	void PostInit(int restart_flag);
	void SetBCTypes();
    void InitFluid();
	void SetBackgroundPressure();
	void InitialProjection();
    void InitialIterations();

    // Member variables for initial conditions
    int probtype = 0;
    Real ic_u = 0.0;
    Real ic_v = 0.0;
    Real ic_w = 0.0;
    Real ic_p = 0.0;

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
    void ComputeDt(int initialisation);
	bool SteadyStateReached();
	void ApplyPredictor();
	void ApplyCorrector();
    void ApplyProjection(Real time, Real scaling_factor);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Convection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

     void incflo_compute_convective_term(Vector< std::unique_ptr<MultiFab> >& conv_u,
                                         Vector< std::unique_ptr<MultiFab> >& conv_r,
                                         Vector< std::unique_ptr<MultiFab> >& conv_t,
                                         Vector< std::unique_ptr<MultiFab> >& vel,
                                         Vector< std::unique_ptr<MultiFab> >& density_in,
                                         Vector< std::unique_ptr<MultiFab> >& tracer_in,
                                         Real time);

    void apply_MAC_projection ( Vector< std::unique_ptr<MultiFab> >& u_mac,
                                Vector< std::unique_ptr<MultiFab> >& v_mac,
                                Vector< std::unique_ptr<MultiFab> >& w_mac,
                                Vector< std::unique_ptr<MultiFab> >& density,
                                Real time);

     void incflo_divergence_plus_redist(const int lev, Vector< std::unique_ptr<MultiFab> >& conv,
                                        Vector< std::unique_ptr<MultiFab> >& fx,
                                        Vector< std::unique_ptr<MultiFab> >& fy,
                                        Vector< std::unique_ptr<MultiFab> >& fz,
                                        int num_comp);

     // Verbosity
     int  mac_mg_verbose = 0;
     int  mac_mg_cg_verbose = 0;

     // Control MLMG behavior
     int  mac_mg_maxiter = 200;
     int  mac_mg_cg_maxiter = 200;
     Real mac_mg_rtol = 1.0e-11;
     Real mac_mg_atol = 1.0e-14;
     std::string mac_bottom_solver_type;

     // Max coarsening level
     int  mac_mg_max_coarsening_level = 100;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Nodal Projection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<NodalProjector> nodal_projector;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<DiffusionOp> diffusion_op;

    static DiffusionType m_diff_type;

    void ComputeDivU(Real time_in);

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

    void MakeEBGeometry();
    bool UpdateEBFactory(int a_lev);

    std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(bool& has_real_walls);

    void make_eb_annulus();
    void make_eb_box();
    void make_eb_cylinder();
    void make_eb_twocylinders();
    void make_eb_regular();
    void make_eb_sphere();
    void make_eb_spherecube();
    void make_eb_tuscan();

    const EB2::Level* eb_level;
    Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;

    // Number of ghost nodes for field arrays.
    const int nghost = 4;

    // These values are required when fluid sees EB -- for now
    const int m_eb_basic_grow_cells = nghost;
    const int m_eb_volume_grow_cells = nghost;
    const int m_eb_full_grow_cells = nghost;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Runtime parameters
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Be verbose?
    int incflo_verbose = 0;

    // Time and time step counters
    Real cur_time = 0.0;
    Real dt = -1.0;
    int nstep = -1;

    // Stop simulation if cur_time reaches stop_time OR nstep reaches max_step
    // OR steady_state = true AND steady_state_tol is reached
    Real stop_time = -1.0;
    int max_step = -1;
    bool steady_state = false;
    Real steady_state_tol = 1.0e-5;

    // Options to control time stepping
    Real cfl = 0.5;
    Real fixed_dt = -1.;

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

    Periodicity p0_periodicity;
    Vector<Real> gp0{Vector<Real>{0.0, 0.0, 0.0}};
    Vector<Real> gravity{Vector<Real>{0.0, 0.0, 0.0}};
    Vector<Real> delp{Vector<Real>{0.0, 0.0, 0.0}};
    Vector<Real> mu_s;
    Real ro_0 = 1.0;

    bool constant_density = true;
    bool advect_tracer  = false;

    // Fluid properties
    std::string fluid_model;
    Real mu = 1.0;
    Real n_0 = 0.0;
    Real tau_0 = 0.0;
    Real papa_reg = 0.0;
    Real eta_0 = 0.0;

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
    Real plot_per = -1.0;
    int last_plt = -1;
    std::string plot_file{"plt"};

    int check_int = -1;
    int last_chk = -1;
    int KE_int = -1;
    std::string check_file{"chk"};
    std::string restart_file{""};

    // Flags for saving fluid data in plot files
    int plt_vel         = 1;
    int plt_gradp       = 1;
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

    // Vector variables to simplify IO -- the vectors are the same for plt/chk
    Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
    Vector<std::string> vecVarsName = {"velx", "vely", "velz", "gpx", "gpy", "gpz"};

    // Scalar variables at cell centers that need to be written to checkfile.
    Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars = {&p, &density, &tracer, &eta};
    Vector<std::string> chkscaVarsName = {"p", "density", "tracer", "eta"};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diagnostics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    Real Norm(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp, int norm_type);
    void PrintMaxValues(Real time);
    void PrintMaxVel(int lev);
    void PrintMaxGp(int lev);
    void CheckForNans(int lev);

    const Real covered_val = 1.0e40;
    const Real boundary_val = 1.0e80;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Data containers
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Boundary conditions
    Vector<std::unique_ptr<IArrayBox>> bc_ilo;
    Vector<std::unique_ptr<IArrayBox>> bc_ihi;
    Vector<std::unique_ptr<IArrayBox>> bc_jlo;
    Vector<std::unique_ptr<IArrayBox>> bc_jhi;
    Vector<std::unique_ptr<IArrayBox>> bc_klo;
    Vector<std::unique_ptr<IArrayBox>> bc_khi;

    // Primary variables
    Vector<std::unique_ptr<MultiFab>> density;
    Vector<std::unique_ptr<MultiFab>> density_o;
    Vector<std::unique_ptr<MultiFab>> tracer;
    Vector<std::unique_ptr<MultiFab>> tracer_o;
    Vector<std::unique_ptr<MultiFab>> vel;
    Vector<std::unique_ptr<MultiFab>> vel_o;
    Vector<std::unique_ptr<MultiFab>> p;
    Vector<std::unique_ptr<MultiFab>> p0;
    Vector<std::unique_ptr<MultiFab>> gp;

    // Derived variables TODO: to save memory, would only need to have 2 temporary variables
    Vector<std::unique_ptr<MultiFab>> eta;
    Vector<std::unique_ptr<MultiFab>> eta_old;
    Vector<std::unique_ptr<MultiFab>> strainrate;
    Vector<std::unique_ptr<MultiFab>> vort;
    Vector<std::unique_ptr<MultiFab>> drag;
    Vector<std::unique_ptr<MultiFab>> divu;

    // Helper variables
    Vector<std::unique_ptr<MultiFab>> conv_u;
    Vector<std::unique_ptr<MultiFab>> conv_u_old;
    Vector<std::unique_ptr<MultiFab>> conv_r;
    Vector<std::unique_ptr<MultiFab>> conv_r_old;
    Vector<std::unique_ptr<MultiFab>> conv_t;
    Vector<std::unique_ptr<MultiFab>> conv_t_old;

    Vector<std::unique_ptr<MultiFab>> divtau;
    Vector<std::unique_ptr<MultiFab>> divtau_old;

    Vector<std::unique_ptr<MultiFab>> laps;
    Vector<std::unique_ptr<MultiFab>> laps_old;

    Vector<std::unique_ptr<MultiFab>> xslopes_u;
    Vector<std::unique_ptr<MultiFab>> yslopes_u;
    Vector<std::unique_ptr<MultiFab>> zslopes_u;

    Vector<std::unique_ptr<MultiFab>> xslopes_r;
    Vector<std::unique_ptr<MultiFab>> yslopes_r;
    Vector<std::unique_ptr<MultiFab>> zslopes_r;

    Vector<std::unique_ptr<MultiFab>> xslopes_t;
    Vector<std::unique_ptr<MultiFab>> yslopes_t;
    Vector<std::unique_ptr<MultiFab>> zslopes_t;

    Vector<std::unique_ptr<MultiFab>> m_u_mac;
    Vector<std::unique_ptr<MultiFab>> m_v_mac;
    Vector<std::unique_ptr<MultiFab>> m_w_mac;
    
    // integer MultiFab specifying if a finer level exists above
    // value is 1 if a finer level is not above and 0 otherwise
    Vector<std::unique_ptr<iMultiFab>> level_mask;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Not yet classified
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void FillPatchVel    (int lev, Real time, MultiFab& mf);
    void FillPatchDensity(int lev, Real time, MultiFab& mf);
    void FillPatchScalar (int lev, Real time, MultiFab& mf);

    void GetDataVel    (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime);
    void GetDataDensity(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime);
    void GetDataScalar (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime);

    void AverageDown();
    void AverageDownTo(int crse_lev);

    void AllocateArrays(int lev);
    void RegridArrays(int lev);
    void MakeBCArrays();

    Vector<Real> t_old;
    Vector<Real> t_new;
};

#endif
