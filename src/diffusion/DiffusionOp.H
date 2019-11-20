#ifndef DIFFUSION_EQUATION_H_
#define DIFFUSION_EQUATION_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#else
#include <AMReX_MLTensorOp.H>
#endif

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

class DiffusionOp
{
public: 
    // Constructor, with initialisation
    DiffusionOp(amrex::AmrCore* _amrcore, 
#ifdef AMREX_USE_EB
                amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* _ebfactory,
#endif
                std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_lo,
                std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_velbc_hi,
                std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_lo,
                std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_scalbc_hi,
                int _nghost,
                int probtype);

    // Destructor
    ~DiffusionOp();

    // Actually do the setup work here
#ifdef AMREX_USE_EB
    void setup(amrex::AmrCore* _amrcore, 
               amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* _ebfactory);
#else
    void setup(amrex::AmrCore* _amrcore); 
#endif

    // Read input from ParmParse (solver settings)
    void readParameters();

    // Set user-supplied solver settings (must be done every time step)
    void setSolverSettings(amrex::MLMG& solver);

    // Tensor viscous solve 
    void diffuse_velocity(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& vel, 
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro, 
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mu, 
                          amrex::Real dt);

    // Scalar diffusion 
    void diffuse_scalar(amrex::Vector<std::unique_ptr<amrex::MultiFab>>& vel, 
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro, 
                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mu,
                        amrex::Real dt);

    // Tensor viscous term 
    void ComputeDivTau(      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& divtau_out,
                       const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& vel_in,
                       const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro_in,
                       const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mu_in);

    // Scalar diffusive term
    void ComputeLapS(      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& laps_out,
                     const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& vel_in,
                     const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ro_in,
                     const amrex::Vector<amrex::Real> mu_s);

private:

    // AmrCore data 
    amrex::AmrCore* amrcore;

#ifdef AMREX_USE_EB
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory;
#endif

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;          
    amrex::Vector<amrex::DistributionMapping> dmap;

    int nghost = 0; 
    int probtype = 0;
    
    // Boundary conditions
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_velbc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_velbc_hi;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_scalbc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_scalbc_hi;

    // Internal data used in the tensor solve
    //
    // MultiLevel EmbeddedBoundary cell-centered tensor solve
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::MLEBTensorOp> vel_matrix;
#else
    std::unique_ptr<amrex::MLTensorOp> vel_matrix;
#endif

    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel_eb;

    // Internal data used in the matrix solve
    //
    // MultiLevel EmbeddedBoundary cell-centered Laplacian: 
    //
    // ( alpha a - beta div ( b grad ) ) phi = rhs
    //
#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::MLEBABecLap>     scal_matrix;
#else
    std::unique_ptr<amrex::MLABecLaplacian> scal_matrix;
#endif

    // DiffusionOp verbosity
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
