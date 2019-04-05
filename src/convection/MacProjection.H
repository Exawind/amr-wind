#ifndef MAC_PROJECTION_H_
#define MAC_PROJECTION_H_

#include <AMReX_AmrCore.H>

#include <constants.H>

class MacProjection
{
public:
	MacProjection(amrex::AmrCore* a_amrcore,
				  int a_nghost,
				  amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* a_ebfactory, 
                  int a_probtype);

	~MacProjection();

	void set_bcs(amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_ilo,
				 amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_ihi,
				 amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_jlo,
				 amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_jhi,
				 amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_klo,
				 amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& a_bc_khi);

	void apply_projection(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& u,
						  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& v,
						  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& w,
						  const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro, 
                          amrex::Real time, int steady_state);

	void update_internals();

	void set_velocity_bcs(int lev,
						  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& u,
						  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& v,
						  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& w, 
                          amrex::Real time);

    void compute_b_coeff(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro, int lev);

private:
	int m_nghost;
    int m_probtype = 0;

	amrex::AmrCore* m_amrcore;
	amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* m_ebfactory;

	// BCs
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_ilo;
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_ihi;
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_jlo;
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_jhi;
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_klo;
	amrex::Vector<std::unique_ptr<amrex::IArrayBox>>* m_bc_khi;

	amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_divu;
	amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_phi;
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> m_b;

	//
	// Stuff for linear solver
	//
	// LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
	std::array<amrex::LinOpBCType, AMREX_SPACEDIM> m_lobc;
	std::array<amrex::LinOpBCType, AMREX_SPACEDIM> m_hibc;

	bool verbose = false;
	int mg_verbose = 0;

	amrex::Real mg_rtol = 1.0e-11;
	amrex::Real mg_atol = 1.0e-14;

    // What solver to use as the bottom solver in the MLMG solves.
    std::string bottom_solver_type;


	void read_inputs();

	//
	// Stuff to compute norms of EB-MultiFabs
	//
	amrex::Real norm0(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mf, int lev);
};

#endif
