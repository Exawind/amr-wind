#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_iMultiFab.H>

#include <eb_if.H>
#include <DiffusionEquation.H>
#include <MacProjection.H>
#include <PoissonEquation.H>


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

private:
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

    // Post-initialization: set BCs, apply ICs, initial velocity projection, pressure iterations
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

	void ComputeUGradU(Vector<std::unique_ptr<MultiFab>>& conv,
					   Vector<std::unique_ptr<MultiFab>>& vel, 
                       Real time);
	void ComputeVelocitySlopes(int lev, MultiFab& Sborder);
	void ComputeVelocityAtFaces(Vector<std::unique_ptr<MultiFab>>& vel, Real time);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Diffusion
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ComputeDivTau(int lev, MultiFab& divtau, Vector<std::unique_ptr<MultiFab>>& vel);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Derived quantities
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void UpdateDerivedQuantities();
	void ComputeDivU(Real time);
	void ComputeStrainrate();
	void ComputeVorticity();
	void ComputeViscosity();
    void ComputeDrag();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Boundary conditions
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void FillScalarBC();
	void FillVelocityBC(Real time, int extrap_dir_bcs);

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

	const EB2::Level* eb_level;
	Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;

	// Number of ghost nodes for field arrays.
    // TODO: can go back to 4 now I think
	const int nghost = 5;

	// These values are required when fluid sees EB -- for now
	const int m_eb_basic_grow_cells = nghost;
	const int m_eb_volume_grow_cells = nghost;
	const int m_eb_full_grow_cells = nghost;

    // Enforce inhomogeneous velocity BC on EB (just cylinder for now)
    Real cyl_speed = 0.0;

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

    // AMR / refinement settings 
	int refine_cutcells = 1;
    int regrid_int = -1;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: Physics
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	Periodicity p0_periodicity;
    Vector<Real> gp0{Vector<Real>{0.0, 0.0, 0.0}};
    Vector<Real> gravity{Vector<Real>{0.0, 0.0, 0.0}};
    Vector<Real> delp{Vector<Real>{0.0, 0.0, 0.0}};
    Real ro_0 = 1.0;

    // Fluid properties
    std::string fluid_model;
    Real mu = 1.0;
    Real n = 0.0;
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
    void WritePlotFile() const;
    void ReadCheckpointFile();

    // Member variables for I/O
    int plot_int = -1;
    Real plot_per = -1.0;
    int last_plt = -1;
    std::string plot_file{"plt"};

    int check_int = -1;
    int last_chk = -1;
    std::string check_file{"chk"};
    std::string restart_file{""};

    // Flags for saving fluid data in plot files
    int plt_vel         = 1;
    int plt_gradp       = 0;
    int plt_rho         = 0;
    int plt_p           = 0;
    int plt_eta         = 1;
    int plt_vort        = 1;
    int plt_strainrate  = 1;
    int plt_stress      = 0;
    int plt_divu        = 0;
    int plt_vfrac       = 1;

    // Total number of variables to write in plot file
    int pltVarCount = 0;



	// Vector variables to simplify IO -- the vectors are the same for plt/chk
	Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
	Vector<std::string> vecVarsName = {"velx", "vely", "velz", "gpx", "gpy", "gpz"};

    // Scalar variables at cell centers that need to be written to checkfile. 
    Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars = {&p, &ro, &eta};
    Vector<std::string> chkscaVarsName = {"p", "ro", "eta"};
    
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

    // Instances of self-written classes for MAC projection, diffusion + Poisson equations 
	std::unique_ptr<MacProjection> mac_projection;
	std::unique_ptr<DiffusionEquation> diffusion_equation;
	std::unique_ptr<PoissonEquation> poisson_equation;

    // Boundary conditions
	Vector<std::unique_ptr<IArrayBox>> bc_ilo;
	Vector<std::unique_ptr<IArrayBox>> bc_ihi;
	Vector<std::unique_ptr<IArrayBox>> bc_jlo;
	Vector<std::unique_ptr<IArrayBox>> bc_jhi;
	Vector<std::unique_ptr<IArrayBox>> bc_klo;
	Vector<std::unique_ptr<IArrayBox>> bc_khi;

    // Primary variables
	Vector<std::unique_ptr<MultiFab>> ro;
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
    Vector<std::unique_ptr<MultiFab>> conv; 
    Vector<std::unique_ptr<MultiFab>> conv_old; 
    Vector<std::unique_ptr<MultiFab>> divtau; 
    Vector<std::unique_ptr<MultiFab>> divtau_old; 
	Vector<std::unique_ptr<MultiFab>> xslopes;
	Vector<std::unique_ptr<MultiFab>> yslopes;
	Vector<std::unique_ptr<MultiFab>> zslopes;
	Vector<std::unique_ptr<MultiFab>> m_u_mac;
	Vector<std::unique_ptr<MultiFab>> m_v_mac;
	Vector<std::unique_ptr<MultiFab>> m_w_mac;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Not yet classified
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void FillPatchVel(int lev, Real time, MultiFab& mf, int icomp, int ncomp);
    void GetDataVel(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime);

	void AverageDown();
	void AverageDownTo(int crse_lev);

	void AllocateArrays(int lev);
	void RegridArrays(int lev);
    void MakeBCArrays();

     Vector<Real> t_old;
     Vector<Real> t_new;
};

#endif
