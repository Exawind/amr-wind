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

class incflo_level : public AmrCore
{
public:
	incflo_level();
	~incflo_level();

	void InitParams();

	void Init(int lev, Real time);

	void InitLevelData(int lev, Real time);

	void PostInit(int lev,
				  Real& dt,
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

	void
	WritePlotFile(std::string& plot_file_name, int nstep = 0, Real dt = 0.0, Real time = 0.0) const;

	void Restart(std::string& restart_chkfile, int* nstep, Real* dt, Real* time, IntVect& Nrep);

	void Regrid(int lev);

	void Advance(
		int lev, int nstep, int steady_state, Real& dt, Real& prev_dt, Real time, Real stop_time);

	void make_mg_bc(int mg_bc[]);

	void make_eb_geometry(int lev);

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_walls(int lev, bool& has_walls);

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> get_real_walls(int lev, bool& has_real_walls);

	std::unique_ptr<EB2::TranslationIF<EB2::PolynomialIF>>
	get_poly(int lev, int max_order, std::string field_prefix);

	std::unique_ptr<EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF>>
	make_wall(int dir, // direction (long edge) of wall
			  Real position,
			  Real height,
			  Real width);

	void WriteEBSurface(int lev);

	void InitIOData();

	static std::string get_load_balance_type();

	void incflo_compute_vort(int lev);

protected:
	//! Tagging cells for refinement
	virtual void
	ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override
	{
		amrex::Abort("ErrorEst: To be implemented");
	}

	//! Make a new level using provided BoxArray and
	//! DistributionMapping and fill with interpolated coarse level
	//! data.  Called by AmrCore::regrid.
	virtual void MakeNewLevelFromCoarse(int lev,
										amrex::Real time,
										const amrex::BoxArray& ba,
										const amrex::DistributionMapping& dm) override
	{
		amrex::Abort("MakeNewLevelFromCoarse: To be implemented");
	}

	//! Remake an existing level using provided BoxArray and
	//! DistributionMapping and fill with existing fine and coarse
	//! data.  Called by AmrCore::regrid.
	virtual void RemakeLevel(int lev,
							 amrex::Real time,
							 const amrex::BoxArray& ba,
							 const amrex::DistributionMapping& dm) override
	{
		amrex::Abort("RemakeLevel: To be implemented");
	}

	//! Delete level data.  Called by AmrCore::regrid.
	virtual void ClearLevel(int lev) override
	{
		amrex::Abort("ClearLevel: To be implemented");
	}

	void incflo_calc_trd_and_tau(int lev);

	void incflo_init_fluid(
		int lev, int is_restarting, Real time, Real& dt, Real stop_time, int steady_state);
	void incflo_set_bc0(int lev);
	void incflo_set_p0(int lev);

	void incflo_project_velocity(int lev);
	void
	incflo_initial_iterations(int lev, amrex::Real dt, amrex::Real stop_time, int steady_state);
	void incflo_apply_projection(int lev, amrex::Real scaling_factor, bool proj_2);
	void incflo_set_scalar_bcs(int lev);
	void incflo_set_velocity_bcs(int lev, int extrap_dir_bcs);

	void AllocateArrays(int lev);
	void RegridArrays(int lev);

	void MakeBCArrays();

private:
	void make_eb_general(int lev);
	void make_eb_annulus(int lev);
	void make_eb_box(int lev);
	void make_eb_cylinder(int lev);
	void make_eb_hopper(int lev);
	void make_eb_hourglass(int lev);
	void make_eb_clr(int lev);
	void make_eb_clr_riser(int lev);
	void make_eb_regular(int lev);

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
	void solve_poisson_equation(int lev,
								Vector<Vector<std::unique_ptr<MultiFab>>>& b,
								Vector<std::unique_ptr<MultiFab>>& phi,
								Vector<std::unique_ptr<MultiFab>>& rhs,
								int bc_lo[],
								int bc_hi[],
								MultiFab& fluxes);

	// Solve implicit diffusion
	void solve_diffusion_equation(int lev,
								  Vector<Vector<std::unique_ptr<MultiFab>>>& b,
								  Vector<std::unique_ptr<MultiFab>>& phi,
								  Vector<std::unique_ptr<MultiFab>>& rhs,
								  int bc_lo[],
								  int bc_hi[],
								  amrex::Real dt);

	void incflo_compute_bcoeff_ppe(int lev);
	void incflo_compute_bcoeff_diff(int lev);

	void check_for_nans(int lev);

	void incflo_compute_dt(
		int lev, amrex::Real time, amrex::Real stop_time, int steady_state, amrex::Real& dt);

	void incflo_compute_ugradu_predictor(int lev,
										 amrex::MultiFab& conv,
										 Vector<std::unique_ptr<MultiFab>>& vel);
	void incflo_compute_ugradu_corrector(int lev,
										 amrex::MultiFab& conv,
										 Vector<std::unique_ptr<MultiFab>>& vel);
	void
	incflo_compute_divtau(int lev, amrex::MultiFab& divtau, Vector<std::unique_ptr<MultiFab>>& vel);

	int steady_state_reached(int lev, Real dt);

	void incflo_apply_predictor(int lev,
								amrex::MultiFab& conv_old,
								amrex::MultiFab& divtau_old,
								amrex::Real dt,
								bool proj_2);
	void incflo_apply_corrector(int lev,
								amrex::MultiFab& conv_old,
								amrex::MultiFab& divtau_old,
								amrex::Real dt,
								bool proj_2);

	void incflo_diffuse_velocity(int lev, amrex::Real dt);

	void incflo_compute_intermediate_velocity(int lev, amrex::Real dt);

	void
	incflo_apply_forcing_terms(int lev, amrex::Real dt, Vector<std::unique_ptr<MultiFab>>& vel);
	void incflo_compute_velocity_slopes(int lev, Vector<std::unique_ptr<MultiFab>>& vel);
	void incflo_compute_velocity_at_faces(int lev, Vector<std::unique_ptr<MultiFab>>& vel);

	void incflo_compute_diveu(int lev);

	void incflo_extrap_pressure(int lev, std::unique_ptr<amrex::MultiFab>& p);

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

   void incflo_average_cc_to_fc ( int lev, const MultiFab& cc,
                                 Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc );
   
	Vector<int> istep = {1};

	// Unit vectors in Cartesian space
	static amrex::IntVect e_x;
	static amrex::IntVect e_y;
	static amrex::IntVect e_z;

	IArrayBox bc_ilo;
	IArrayBox bc_ihi;
	IArrayBox bc_jlo;
	IArrayBox bc_jhi;
	IArrayBox bc_klo;
	IArrayBox bc_khi;

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
	Vector<std::unique_ptr<MultiFab>> mu;

	// Cell-based
	Vector<std::unique_ptr<MultiFab>> lambda;
	Vector<std::unique_ptr<MultiFab>> trD;
	Vector<std::unique_ptr<MultiFab>> vort;

	// These are multi-component multifabs
	Vector<std::unique_ptr<MultiFab>> xslopes;
	Vector<std::unique_ptr<MultiFab>> yslopes;
	Vector<std::unique_ptr<MultiFab>> zslopes;

	// div (u)
	Vector<std::unique_ptr<MultiFab>> diveu;

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
	// Vector< std::unique_ptr<MultiFab> > m_uo_mac;
	// Vector< std::unique_ptr<MultiFab> > m_vo_mac;
	// Vector< std::unique_ptr<MultiFab> > m_wo_mac;

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
    Real mu_0 = 1.0;

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

	int nodal_pressure = 0;
	int explicit_diffusion = 1;

	int ugradu_type = 3;

	// Variables to simplify IO -- the vectors are the same for plt/chk
	Vector<const Vector<std::unique_ptr<MultiFab>>*> vectorVars;
	Vector<std::string> vecVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> chkscalarVars;
	Vector<std::string> pltscaVarsName;

	Vector<const Vector<std::unique_ptr<MultiFab>>*> pltscalarVars;
	Vector<std::string> chkscaVarsName;

	/****************************************************************************
     *                                                                          *
     * EB Data (including level-set data)                                       *
     *                                                                          *
     ****************************************************************************/

	const EB2::Level* eb_level_fluid;

	Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>> ebfactory;
	std::unique_ptr<MultiFab> dummy;

    void incflo_update_ebfactory(int a_lev);
};

inline std::string incflo_level::get_load_balance_type()
{
	return load_balance_type;
}

#endif
