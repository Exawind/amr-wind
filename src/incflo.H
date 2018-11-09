#ifndef INCFLO_H_
#define INCFLO_H_

#include <AMReX_iMultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>

#include <eb_if.H>
#include <MacProjection.H>
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
    
    void ReadParameters();

	void InitParams();

	void Init();

	void InitLevelData(Real time);

	void PostInit(Real& dt,
				  Real time,
				  int nstep,
				  int restart_flag,
				  Real stop_time,
				  int steady_state);

	void ResizeArrays();

    void incflo_init_fluid(int is_restarting, 
                           Real time, 
                           Real& dt, 
                           Real stop_time, 
                           int steady_state);

	void incflo_set_bc0();
	void incflo_set_p0();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Input / Output
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void WriteCheckPointFile(std::string& check_file_name,
							 int nstep = 0,
							 Real dt = 0.0,
							 Real time = 0.0) const;

    void WritePlotFile(std::string& plot_file_name, 
                       int nstep = 0, 
                       Real dt = 0.0, 
                       Real time = 0.0) const;

    void Restart(std::string& restart_chkfile, int* nstep, Real* dt, Real* time);

    void WriteHeader(const std::string& name, int nstep, Real dt, Real time, 
                     bool is_checkpoint) const;
	void WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const;
	void WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const;
	void WriteJobInfo(const std::string& dir) const;


    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Simulate!
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void Advance(int nstep, 
                 int steady_state, 
                 Real& dt, 
                 Real& prev_dt, 
                 Real time, 
                 Real stop_time);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Embedded Boundaries
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void make_eb_geometry();

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls(bool& has_walls);

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(bool& has_real_walls);

    std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>> get_poly(int max_order, 
                                                                    std::string field_prefix);

	std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>
	make_wall(int dir, // direction (long edge) of wall
			  Real position,
			  Real height,
			  Real width);

	void WriteEBSurface();

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

	void incflo_compute_strainrate();
	void incflo_compute_vort();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Rheology
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void incflo_compute_viscosity();

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Projection
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

    void incflo_apply_projection(Real time, Real scaling_factor, bool proj_2);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Boundary conditions
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void incflo_set_scalar_bcs();
	void incflo_set_velocity_bcs(Real time, int extrap_dir_bcs);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Not yet classified
    //
    //////////////////////////////////////////////////////////////////////////////////////////////

	void make_mg_bc(int mg_bc[]);

	void InitIOData();

	void incflo_initial_projection();
    void incflo_initial_iterations(Real dt, Real stop_time, int steady_state);

	void FillPatch(int lev, MultiFab& mf, MultiFab& cmf, MultiFab &fmf, int icomp, int ncomp);

	void AllocateArrays(int lev);
	void RegridArrays(int lev);

    void MakeBCArrays();

	void ReMakeNewLevelFromScratch(int lev,
								   const BoxArray& new_grids,
								   const DistributionMapping& new_dmap);

	static void GotoNextLine(std::istream& is);

	// Projection-related methods
	void solve_poisson_equation(Vector<Vector<std::unique_ptr<MultiFab>>>& b,
								Vector<std::unique_ptr<MultiFab>>& phi,
								Vector<std::unique_ptr<MultiFab>>& rhs,
								Vector<std::unique_ptr<MultiFab>>& fluxes,
                                int bc_lo[], int bc_hi[]);

	// Solve implicit diffusion
	void solve_diffusion_equation(Vector<Vector<std::unique_ptr<MultiFab>>>& b,
								  Vector<std::unique_ptr<MultiFab>>& phi,
								  Vector<std::unique_ptr<MultiFab>>& rhs,
                                  int bc_lo[], int bc_hi[],
								  Real dt);

	void incflo_compute_bcoeff_ppe();
	void incflo_compute_bcoeff_diff();

	void check_for_nans(int lev);

	void incflo_compute_dt(Real time, 
                           Real stop_time, 
                           int steady_state, 
                           int initialisation, 
                           Real& dt);

	void incflo_compute_ugradu_predictor(Vector<std::unique_ptr<MultiFab>>& conv,
										 Vector<std::unique_ptr<MultiFab>>& vel, 
                                         Real time);
	void incflo_compute_ugradu_corrector(Vector<std::unique_ptr<MultiFab>>& conv,
										 Vector<std::unique_ptr<MultiFab>>& vel, 
                                         Real time);
    void incflo_compute_divtau(int lev, 
                               MultiFab& divtau, 
                               Vector<std::unique_ptr<MultiFab>>& vel);

	int steady_state_reached(Real dt, int iter);

	void incflo_apply_predictor(Vector<std::unique_ptr<MultiFab>>& conv_old,
								Vector<std::unique_ptr<MultiFab>>& divtau_old,
                                Real time, Real dt, bool proj_2);
	void incflo_apply_corrector(Vector<std::unique_ptr<MultiFab>>& conv_old,
								Vector<std::unique_ptr<MultiFab>>& divtau_old,
                                Real time, Real dt, bool proj_2);

	void incflo_diffuse_velocity(Real time, Real dt);

	void incflo_compute_velocity_slopes(int lev, Vector<std::unique_ptr<MultiFab>>& vel);
	void incflo_compute_velocity_at_faces(int lev, Vector<std::unique_ptr<MultiFab>>& vel);

	void incflo_compute_divu(Real time);

	void incflo_extrap_pressure(int lev, std::unique_ptr<MultiFab>& p);

	void incflo_print_max_vel(int lev);

	// Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
	BoxArray MakeBaseGrids() const;
	void ChopGrids(const Box& domain, BoxArray& ba, int target_size) const;

	void incflo_set_bc_type(int lev);

	void fill_mf_bc(int lev, MultiFab& mf);

    Real incflo_norm1(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp);
    Real incflo_norm1(MultiFab& mf, int lev, int comp);
    Real incflo_norm0(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp);
	Real incflo_norm0(MultiFab& mf, int lev, int comp);

    void incflo_average_cc_to_fc (int lev, 
                                  const MultiFab& cc,
                                  Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc);

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
    // but will never change unless we use subcycling. 
    Vector<int> istep = {1};

	Vector<std::unique_ptr<IArrayBox>> bc_ilo;
	Vector<std::unique_ptr<IArrayBox>> bc_ihi;
	Vector<std::unique_ptr<IArrayBox>> bc_jlo;
	Vector<std::unique_ptr<IArrayBox>> bc_jhi;
	Vector<std::unique_ptr<IArrayBox>> bc_klo;
	Vector<std::unique_ptr<IArrayBox>> bc_khi;

	Vector<std::unique_ptr<iMultiFab>> flag;

	// Gas pressure fraction
	Vector<std::unique_ptr<MultiFab>> p;

	// Gas density
	Vector<std::unique_ptr<MultiFab>> ro;

	// Gas velocity
	Vector<std::unique_ptr<MultiFab>> vel;
	Vector<std::unique_ptr<MultiFab>> vel_o;

	// Base state pressure
	Vector<std::unique_ptr<MultiFab>> p0;

	// Pressure gradients
	Vector<std::unique_ptr<MultiFab>> gp;
	Vector<std::unique_ptr<MultiFab>> gp0;

	// Molecular viscosity
	Vector<std::unique_ptr<MultiFab>> eta;

	// Strain-rate magnitude
	Vector<std::unique_ptr<MultiFab>> strainrate;

    // Vorticity
	Vector<std::unique_ptr<MultiFab>> vort;

	// These are multi-component multifabs
	Vector<std::unique_ptr<MultiFab>> xslopes;
	Vector<std::unique_ptr<MultiFab>> yslopes;
	Vector<std::unique_ptr<MultiFab>> zslopes;

	// div (u)
	Vector<std::unique_ptr<MultiFab>> divu;

	// RHS for implicit diffusion solve
	Vector<std::unique_ptr<MultiFab>> rhs_diff;

	// Solution for implicit diffusion solve
	Vector<std::unique_ptr<MultiFab>> phi_diff;

	//
	Vector<std::unique_ptr<MultiFab>> fp;

	// Coefficients b of div(b*grad(p))
	Vector<Vector<std::unique_ptr<MultiFab>>> bcoeff;

	// Coefficients b in implicit diffusion solve
	Vector<Vector<std::unique_ptr<MultiFab>>> bcoeff_diff;

	// Pressure increment
	Vector<std::unique_ptr<MultiFab>> phi;

	// MAC velocity
	Vector<std::unique_ptr<MultiFab>> m_u_mac;
	Vector<std::unique_ptr<MultiFab>> m_v_mac;
	Vector<std::unique_ptr<MultiFab>> m_w_mac;

	// MAC projection object
	std::unique_ptr<MacProjection> mac_projection;

    // What solver to use as the bottom solver in the MLMG solves.
    std::string bottom_solver_type;

    // Time and time step counters
    // TODO: time should be called something else (like t)
    // TODO: remove unnecessary function parameters since there are now member functions
	Real time = 0.0L;
	int nstep = 0;
	Real dt = -1;

    // Be verbose? 
    int verbose = 0;

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

	// Options to control MLMG behavior
    // TODO: get this out of here , should not be member functions
	int mg_verbose = 0;
	int mg_cg_verbose = 0;
	int mg_max_iter = 100;
	int mg_cg_maxiter = 100;
	int mg_max_fmg_iter = 0;
	Real mg_rtol = 1.0e-11;
	Real mg_atol = 1.0e-14;

	// Tolerance to check for steady state
	Real steady_state_tol = 1.0e-5;

	// Options to control time stepping
	Real cfl = 0.5;
	Real fixed_dt = -1.;

	Periodicity p0_periodicity;

	//
	// Here we set the number of ghost nodes of the field arrays.
	// Since this is  a class attribute, it will be propagated to
	// all the member functions.
	//
	// DO NOT USE  A LOCAL VARIABLE WITH THE SAME NAME
	//

	// These values are required when fluid sees EB -- for now
	const int nghost = 5;

	const int m_eb_basic_grow_cells = nghost;
	const int m_eb_volume_grow_cells = nghost;
	const int m_eb_full_grow_cells = nghost;

    int nodal_pressure = 1;
	int explicit_diffusion = 1;
	int refine_cutcells = 1;

	// Variables to simplify IO -- the vectors are the same for plt/chk
	Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
	Vector<std::string> vecVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars;
	Vector<std::string> pltscaVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> pltscalarVars;
	Vector<std::string> chkscaVarsName;

    // Declare and initialise variables
    Real stop_time = -1.0;
    int max_step = -1;
    bool steady_state = false;

    int check_int = -1;
    int last_chk = -1;
    std::string check_file{"chk"};
    std::string restart_file{""};

    int plot_int = -1;
    int last_plt = -1;
    std::string plot_file{"plt"};
    bool write_eb_surface = false;

    int regrid_int = -1;

	/****************************************************************************
     *                                                                          *
     * EB Data                                                                  *
     *                                                                          *
     ****************************************************************************/

	const EB2::Level* eb_level_fluid;

	Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;

    void incflo_update_ebfactory(int a_lev);
};

#endif
