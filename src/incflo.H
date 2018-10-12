#ifndef INCFLO_LEVEL_H_
#define INCFLO_LEVEL_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB2_IF_Translation.H>

#include <eb_if.H>
#include <MacProjection.H>

class incflo : public AmrCore
{
public:
	incflo();
	~incflo();

	void InitParams();

	void Init(Real time);

	void InitLevelData(Real time);

	void PostInit(Real& dt,
				  Real time,
				  int nstep,
				  int restart_flag,
				  Real stop_time,
				  int steady_state);

	void ResizeArrays();

	void WriteCheckPointFile(std::string& check_file_name,
							 int nstep = 0,
							 Real dt = 0.0,
							 Real time = 0.0) const;

    void WritePlotFile(std::string& plot_file_name, 
                       int nstep = 0, 
                       Real dt = 0.0, 
                       Real time = 0.0) const;

	void Restart(std::string& restart_chkfile, 
                 int* nstep, 
                 Real* dt, 
                 Real* time, 
                 IntVect& Nrep);

	void Regrid();

    void Advance(int nstep, 
                 int steady_state, 
                 Real& dt, 
                 Real& prev_dt, 
                 Real time, 
                 Real stop_time);

	void make_mg_bc(int mg_bc[]);

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

	void InitIOData();

	static std::string get_load_balance_type();

	void incflo_compute_strainrate();
	void incflo_compute_vort();

	void incflo_compute_viscosity();

protected:
	//! Tagging cells for refinement
	virtual void
	ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override
	{
		Abort("ErrorEst: To be implemented");
	}

	//! Make a new level using provided BoxArray and
	//! DistributionMapping and fill with interpolated coarse level
	//! data.  Called by AmrCore::regrid.
	virtual void MakeNewLevelFromCoarse(int lev,
										Real time,
										const BoxArray& ba,
										const DistributionMapping& dm) override
	{
		Abort("MakeNewLevelFromCoarse: To be implemented");
	}

	//! Remake an existing level using provided BoxArray and
	//! DistributionMapping and fill with existing fine and coarse
	//! data.  Called by AmrCore::regrid.
	virtual void RemakeLevel(int lev,
							 Real time,
							 const BoxArray& ba,
							 const DistributionMapping& dm) override
	{
		Abort("RemakeLevel: To be implemented");
	}

	//! Delete level data.  Called by AmrCore::regrid.
	virtual void ClearLevel(int lev) override
	{
		Abort("ClearLevel: To be implemented");
	}

    void incflo_init_fluid(int is_restarting, 
                           Real time, 
                           Real& dt, 
                           Real stop_time, 
                           int steady_state);

	void incflo_set_bc0();
	void incflo_set_p0();

	void incflo_initial_projection();
    void incflo_initial_iterations(Real dt, Real stop_time, int steady_state);

    void incflo_apply_projection(Real dt, Real scaling_factor, bool proj_2);

	void incflo_set_scalar_bcs();
	void incflo_set_velocity_bcs(Real time, int extrap_dir_bcs);

	void FillPatch(int lev, MultiFab& mf, MultiFab& cmf, MultiFab &fmf, int icomp, int ncomp);

	void AllocateArrays(int lev);
	void RegridArrays(int lev);

private:
	void make_eb_general();
	void make_eb_annulus();
	void make_eb_box();
	void make_eb_cylinder();
	void make_eb_hopper();
	void make_eb_hourglass();
	void make_eb_clr();
	void make_eb_clr_riser();
	void make_eb_regular();

    void MakeBCArrays();

	void MakeNewLevelFromScratch(int lev,
								 Real time,
								 const BoxArray& new_grids,
								 const DistributionMapping& new_dmap) override;

	void ReMakeNewLevelFromScratch(int lev,
								   const BoxArray& new_grids,
								   const DistributionMapping& new_dmap);

	void
	WriteHeader(const std::string& name, int nstep, Real dt, Real time, bool is_checkpoint) const;
	void WritePlotHeader(const std::string& name, int nstep, Real dt, Real time) const;
	void WriteCheckHeader(const std::string& name, int nstep, Real dt, Real time) const;

	void WriteJobInfo(const std::string& dir) const;

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
#if 0
    void WritePlotFile() const;

      // which step?
    Vector<int> nsubsteps;  // how many substeps on each level?
#endif

    Real incflo_norm1(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp);
    Real incflo_norm1(MultiFab& mf, int lev, int comp);
    Real incflo_norm0(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp);
	Real incflo_norm0(MultiFab& mf, int lev, int comp);

    void incflo_average_cc_to_fc (int lev, 
                                  const MultiFab& cc,
                                  Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc);
   
    static int nlev; 

	Vector<int> istep = {1};

	// Unit vectors in Cartesian space
	static IntVect e_x;
	static IntVect e_y;
	static IntVect e_z;

	Vector<std::unique_ptr<IArrayBox>> bc_ilo;
	Vector<std::unique_ptr<IArrayBox>> bc_ihi;
	Vector<std::unique_ptr<IArrayBox>> bc_jlo;
	Vector<std::unique_ptr<IArrayBox>> bc_jhi;
	Vector<std::unique_ptr<IArrayBox>> bc_klo;
	Vector<std::unique_ptr<IArrayBox>> bc_khi;

	Vector<std::unique_ptr<iMultiFab>> flag;

	// Gas pressure fraction
	Vector<std::unique_ptr<MultiFab>> p;
	Vector<std::unique_ptr<MultiFab>> p_o;

	// Gas density
	Vector<std::unique_ptr<MultiFab>> ro;
	Vector<std::unique_ptr<MultiFab>> ro_o;

	// Gas velocity
	Vector<std::unique_ptr<MultiFab>> vel;
	Vector<std::unique_ptr<MultiFab>> vel_o;

	// Base state pressure
	Vector<std::unique_ptr<MultiFab>> p0;

	// Pressure correction
	Vector<std::unique_ptr<MultiFab>> pp;

	// Pressure gradients
	Vector<std::unique_ptr<MultiFab>> gp;
	Vector<std::unique_ptr<MultiFab>> gp0;

	// Molecular viscosity
	Vector<std::unique_ptr<MultiFab>> eta;

	// Cell-based
	Vector<std::unique_ptr<MultiFab>> strainrate;
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

	// used if load_balance_type == "KnapSack"
	Vector<std::unique_ptr<MultiFab>> fluid_cost;

    // Be verbose? 
    int verbose = 0;

    // Physics 
    static Vector<Real> gravity;
    Real ro_0 = 1.0;

    // Fluid properties
    std::string fluid_model;
    Real mu = 1.0;
    Real n = 0.0;
    Real tau_0 = 0.0;
    Real papa_reg = 0.0;
    Real eta_0 = 0.0;

	// Options to control MLMG behavior
	int mg_verbose = 0;
	int mg_cg_verbose = 0;
	int mg_max_iter = 100;
	int mg_cg_maxiter = 100;
	int mg_max_fmg_iter = 0;
	Real mg_rtol = 1.0e-11;
	Real mg_atol = 1.0e-14;

	// Tolerance to check for steady state
	Real steady_state_tol = 1.0e-5;

	static std::string load_balance_type;
	static std::string knapsack_weight_type;

	// Options to control time stepping
	Real cfl = 0.5;
	Real fixed_dt = -1.;

	static EBSupport m_eb_support_level;

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

	int explicit_diffusion = 1;

	// Variables to simplify IO -- the vectors are the same for plt/chk
	Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
	Vector<std::string> vecVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars;
	Vector<std::string> pltscaVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> pltscalarVars;
	Vector<std::string> chkscaVarsName;

	/****************************************************************************
     *                                                                          *
     * EB Data                                                                  *
     *                                                                          *
     ****************************************************************************/

	const EB2::Level* eb_level_fluid;

	Vector<std::unique_ptr<EBFArrayBoxFactory>> ebfactory;

    void incflo_update_ebfactory(int a_lev);
};

inline std::string incflo::get_load_balance_type()
{
	return load_balance_type;
}

#endif
