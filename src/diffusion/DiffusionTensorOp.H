#ifndef DIFFUSION_TENSOR_OP_H_
#define DIFFUSION_TENSOR_OP_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBTensorOp.H>
#endif
#include <AMReX_MLTensorOp.H>

//
// Solver for the implicit part of the diffusion equation: 
//
// ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
//
//  rhs = u + dt * ( conv + div ( eta (grad u)^T ) / rho + g - grad(p + p0) / rho )
//

class incflo;

class DiffusionTensorOp
{
public:
    DiffusionTensorOp (incflo* a_incflo);

    void diffuse_velocity (amrex::Vector<amrex::MultiFab*> const& velocity,
                           amrex::Vector<amrex::MultiFab*> const& density,
                           amrex::Vector<amrex::MultiFab const*> const& eta,
                           amrex::Real t, amrex::Real dt);

    void compute_divtau (amrex::Vector<amrex::MultiFab*> const& divtau,
                         amrex::Vector<amrex::MultiFab const*> const& velocity,
                         amrex::Vector<amrex::MultiFab const*> const& density,
                         amrex::Vector<amrex::MultiFab const*> const& eta,
                         amrex::Real t);

private:

    void readParameters ();

    incflo* m_incflo;

#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::MLEBTensorOp> m_eb_solve_op;
    std::unique_ptr<amrex::MLEBTensorOp> m_eb_apply_op;
#endif
    std::unique_ptr<amrex::MLTensorOp> m_reg_solve_op;
    std::unique_ptr<amrex::MLTensorOp> m_reg_apply_op;

    // DiffusionOp verbosity
    int m_verbose = 0;

    // Options to control MLMG behavior
    int m_mg_verbose = 0;
    int m_mg_cg_verbose = 0;
    int m_mg_max_iter = 100;
    int m_mg_cg_maxiter = 100;
    int m_mg_max_fmg_iter = 0;
    int m_mg_max_coarsening_level = 100;
    int m_mg_maxorder = 2;
    amrex::Real m_mg_rtol = 1.0e-11;
    amrex::Real m_mg_atol = 1.0e-14;
    std::string m_bottom_solver_type = "bicgstab";
};

#endif
