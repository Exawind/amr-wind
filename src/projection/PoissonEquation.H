#ifndef POISSON_EQUATION_H_
#define POISSON_EQUATION_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

// TODO: DOCUMENTATION

class PoissonEquation
{
public: 
    // Constructor, with initialisation
    PoissonEquation(amrex::AmrCore* _amrcore, 
                    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* _ebfactory,
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_ilo, 
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_ihi, 
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_jlo, 
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_jhi, 
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_klo, 
                    amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_khi,
                    int _nghost);

    // Destructor
    ~PoissonEquation();

    // Read input from ParmParse (solver settings)
    void readParameters();

    // Update internals if AmrCore changes (e.g. after regrid)
    void updateInternals(amrex::AmrCore* amrcore_in, 
                         amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory_in);

    // Set user-supplied solver settings (must be done every time step)
    void setSolverSettings(amrex::MLMG& solver);

    // Solve the Poisson equation, put results in phi and fluxes
    void solve(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& phi, 
               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& fluxes,
               const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro,
               const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& divu);

private:
    // AmrCore data 
    amrex::AmrCore* amrcore;
	amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory;
    int nghost; 

    // Internal data used in the matrix solve
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> sigma;
    amrex::MLNodeLaplacian matrix;

    // Boundary conditions
    int bc_lo[3], bc_hi[3];

    // PoissonEquation verbosity
	int verbose = 0;

	// Options to control solver behavior
	int mg_verbose = 0;
	int mg_cg_verbose = 0;
	int mg_max_iter = 100;
	int mg_cg_maxiter = 100;
	int mg_max_fmg_iter = 0;
    int mg_max_coarsening_level = 100;
    amrex::Real mg_rtol = 1.0e-11;
    amrex::Real mg_atol = 1.0e-14;
    std::string bottom_solver_type = "bicgcg";
};


#endif
