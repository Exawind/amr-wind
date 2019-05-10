#ifndef DIFFUSION_EQUATION_H_
#define DIFFUSION_EQUATION_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>

//
// Solver for the implicit part of the diffusion equation: 
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
// Input "vel" holds the right hand side of this equation:
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//
// This is replaced by the solution of the matrix solve. 
// Output "vel" therefore holds u*
//
// Note: we actually solve the above equation multiplied by the density ro.
//

class DiffusionEquation
{
public: 
    // Constructor, with initialisation
    DiffusionEquation(amrex::AmrCore* _amrcore, 
                      amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* _ebfactory,
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_ilo, 
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_ihi, 
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_jlo, 
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_jhi, 
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_klo, 
                      amrex::Vector<std::unique_ptr<amrex::IArrayBox>>& bc_khi,
                      int _nghost, amrex::Real cyl_speed);

    // Destructor
    ~DiffusionEquation();

    // Read input from ParmParse (solver settings)
    void readParameters();

    // Update internals if AmrCore changes (e.g. after regrid)
    void updateInternals(amrex::AmrCore* amrcore_in, 
                         amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory_in);

    // Set user-supplied solver settings (must be done every time step)
    void setSolverSettings(amrex::MLMG& solver);

    // Solve the diffusion equation, update vel
    void solve(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& vel, 
               const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro, 
               const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& eta, 
               amrex::Real dt);

private:
    // AmrCore data 
    amrex::AmrCore* amrcore;
	amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory;
    int nghost; 
    amrex::Real cyl_speed = 0.0;

    // Internal data used in the matrix solve
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian: 
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
    amrex::MLEBABecLap matrix;
    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> ueb;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> veb;

    // Boundary conditions
    int bc_lo[3], bc_hi[3];

    // DiffusionEquation verbosity
	int verbose = 0;

	// Options to control MLMG behavior
	int mg_verbose = 0;
	int mg_cg_verbose = 0;
	int mg_max_iter = 100;
	int mg_cg_maxiter = 100;
	int mg_max_fmg_iter = 0;
	int mg_max_coarsening_level = 100;
    amrex::Real mg_rtol = 1.0e-11;
    amrex::Real mg_atol = 1.0e-14;
    std::string bottom_solver_type = "bicgstab";
};


#endif
