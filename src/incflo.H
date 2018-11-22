#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>

#include <eb_if.H>
#include <DiffusionEquation.H>
#include <MacProjection.H>
#include <PoissonEquation.H>
// Note: a lot of stuff is included through MacProjection.H...


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

private:
    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Initialization
    //
    //////////////////////////////////////////////////////////////////////////////////////////////
    
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

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Simulate!
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void Advance();
    void ComputeDt(int initialisation);
	bool SteadyStateReached();
	void ApplyPredictor(Vector<std::unique_ptr<MultiFab>>& conv_old,
						Vector<std::unique_ptr<MultiFab>>& divtau_old);
	void ApplyCorrector(Vector<std::unique_ptr<MultiFab>>& conv_old,
						Vector<std::unique_ptr<MultiFab>>& divtau_old);


    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Embedded Boundaries
    // TODO: How much of this really needs to be member functions of incflo? 
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void MakeEBGeometry();
    void UpdateEBFactory(int a_lev);
	void WriteEBSurface();

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls(bool& has_walls);
	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(bool& has_real_walls);
    std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>> get_poly(int max_order, std::string field_prefix);
	std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>> make_wall(int dir, Real position, Real height, Real width);

	void make_eb_general();
	void make_eb_annulus();
	void make_eb_box();
	void make_eb_cylinder();
	void make_eb_hopper();
	void make_eb_hourglass();
	void make_eb_clr();
	void make_eb_clr_riser();
	void make_eb_regular();

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

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Input / Output
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Write header, which is the same for plot and checkpoint files 
    void WriteHeader(const std::string& name, bool is_checkpoint) const;

    // Write job info (compilers, #processors, git hash etc etc)
	void WriteJobInfo(const std::string& dir) const;

    // Write checkpoint file to restart from 
    void WriteCheckPointFile() const;

    // Write plot file for visualisation of results
    void WritePlotFile() const;

    // Read checkpoint file and restart simulation from it
    void ReadCheckpointFile();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Projection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void ApplyProjection(Real time, Real scaling_factor);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Boundary conditions
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void FillScalarBC(int lev, MultiFab& mf);
	void FillVelocityBC(Real time, int extrap_dir_bcs);

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

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Not yet classified
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void make_mg_bc(int mg_bc[]);

	void FillPatch(int lev, Real time, MultiFab& mf, MultiFab& cmf, MultiFab &fmf, int icomp, int ncomp);

	void AverageDown();
	void AverageDownTo(int crse_lev);

	void AllocateArrays(int lev);
	void RegridArrays(int lev);
    void MakeBCArrays();


    // DIFFUSION
    void ComputeDivTau(int lev, MultiFab& divtau, Vector<std::unique_ptr<MultiFab>>& vel);
	void DiffuseVelocity(Real time);

    // CONVECTION
	void ComputeUGradU(Vector<std::unique_ptr<MultiFab>>& conv,
					   Vector<std::unique_ptr<MultiFab>>& vel, 
                       Real time);
	void ComputeVelocitySlopes(int lev, Vector<std::unique_ptr<MultiFab>>& vel);
	void ComputeVelocityAtFaces(int lev, Vector<std::unique_ptr<MultiFab>>& vel);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: MultiFabs and IArrayBoxes existing on each level
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	Vector<std::unique_ptr<IArrayBox>> bc_ilo;
	Vector<std::unique_ptr<IArrayBox>> bc_ihi;
	Vector<std::unique_ptr<IArrayBox>> bc_jlo;
	Vector<std::unique_ptr<IArrayBox>> bc_jhi;
	Vector<std::unique_ptr<IArrayBox>> bc_klo;
	Vector<std::unique_ptr<IArrayBox>> bc_khi;
	Vector<std::unique_ptr<iMultiFab>> flag;
	Vector<std::unique_ptr<MultiFab>> p;
	Vector<std::unique_ptr<MultiFab>> ro;
	Vector<std::unique_ptr<MultiFab>> vel;
	Vector<std::unique_ptr<MultiFab>> vel_o;
	Vector<std::unique_ptr<MultiFab>> p0;
	Vector<std::unique_ptr<MultiFab>> gp;
	Vector<std::unique_ptr<MultiFab>> gp0;
	Vector<std::unique_ptr<MultiFab>> eta;
	Vector<std::unique_ptr<MultiFab>> strainrate;
	Vector<std::unique_ptr<MultiFab>> vort;
	Vector<std::unique_ptr<MultiFab>> xslopes;
	Vector<std::unique_ptr<MultiFab>> yslopes;
	Vector<std::unique_ptr<MultiFab>> zslopes;
	Vector<std::unique_ptr<MultiFab>> divu;
	Vector<std::unique_ptr<MultiFab>> m_u_mac;
	Vector<std::unique_ptr<MultiFab>> m_v_mac;
	Vector<std::unique_ptr<MultiFab>> m_w_mac;

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: I/O
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    int plot_int = -1;
    int last_plt = -1;
    std::string plot_file{"plt"};

    int check_int = -1;
    int last_chk = -1;
    std::string check_file{"chk"};
    std::string restart_file{""};

	// Vector variables to simplify IO -- the vectors are the same for plt/chk
	Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
	Vector<std::string> vecVarsName = {"velx", "vely", "velz", "gpx", "gpy", "gpz"};

    // Scalar variables at cell centers that need to be written to plotfile/checkfile. 
    // "volfrac" MUST always be last without any mf associated to it!
    Vector<const Vector<std::unique_ptr<MultiFab>>*> pltscalarVars = {&p, &ro, &eta, &strainrate, &strainrate, &vort, &divu};
    Vector<std::string> pltscaVarsName = {"p", "ro", "eta", "strainrate", "stress", "vort", "divu", "volfrac"};
	Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars = {&p, &ro, &eta};
	Vector<std::string> chkscaVarsName = {"p", "ro", "eta"};

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Member variables: unassigned
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	std::unique_ptr<MacProjection> mac_projection;
	std::unique_ptr<DiffusionEquation> diffusion_equation;
	std::unique_ptr<PoissonEquation> poisson_equation;

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

    // Be verbose? 
    int incflo_verbose = 0;

    // Physics 
    Vector<Real> gravity{Vector<Real>{0.0, 0.0, 0.0}};
    Real ro_0 = 1.0;

    // Fluid properties
    std::string fluid_model;
    Real mu = 1.0;
    Real n = 0.0;
    Real tau_0 = 0.0;
    Real papa_reg = 0.0;
    Real eta_0 = 0.0;

	// Options to control time stepping
	Real cfl = 0.5;
	Real fixed_dt = -1.;

	Periodicity p0_periodicity;

	// Number of ghost nodes for field arrays.
	const int nghost = 5;

	// These values are required when fluid sees EB -- for now
	const int m_eb_basic_grow_cells = nghost;
	const int m_eb_volume_grow_cells = nghost;
	const int m_eb_full_grow_cells = nghost;

	int explicit_diffusion = 1;
	int refine_cutcells = 1;
    int regrid_int = -1;


	const EB2::Level* eb_level_fluid;
	Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;
};

#endif
