#ifndef NODAL_PROJECTION_H
#define NODAL_PROJECTION_H

#include <AMReX.H>
#include <incflo.H>

class NodalProjection
{

public:

    // Constructor
    NodalProjection (const incflo* a_incflo,
                     Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                     std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_lo,
                     std::array<amrex::LinOpBCType,AMREX_SPACEDIM> a_bc_hi);

    // Destructor
    NodalProjection () {};

    void setup (const incflo* a_incflo,
                Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory);

    void project (       Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                   const Vector< std::unique_ptr< amrex::MultiFab > >& a_ro,
                         Real a_time, Real a_scale_factor = 1.0 );

    Vector< const amrex::MultiFab* > getPhi () const {return GetVecOfConstPtrs(m_phi);};

    Vector< const amrex::MultiFab* > getGradPhi () const {return GetVecOfConstPtrs(m_fluxes);};


private:

    bool m_ok = false;

    // Verbosity
    int  m_verbose = 0;
    int  m_mg_verbose = 0;
    int  m_mg_cg_verbose = 0;

    // Control MLMG behavior
    int  m_mg_maxiter = 100;
    int  m_mg_cg_maxiter = 100;
    Real m_mg_rtol = 1.0e-11;
    Real m_mg_atol = 1.0e-14;
    std::string m_bottom_solver_type = "bicgcg";

    // Max coarsening level
    int  m_mg_max_coarsening_level = 100;

    // We need ebfactory, amrcore and some functions defined in incflo.
    // Thus we make this class a friend of incflo' and use its methods.
    // In the future we will decouple those method from inclfo and just pass
    // amrcore and ebfactory to the constructor
    const incflo* m_incflo;

    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>* ebfactory;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    // Boundary conditions
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_lo;
    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>  m_bc_hi;

    // Cell-centered data
    Vector< std::unique_ptr< amrex::MultiFab > >  m_fluxes;
    Vector< std::unique_ptr< amrex::MultiFab > >  m_sigma;

    // Node-centered data
    Vector< std::unique_ptr< amrex::MultiFab > >  m_phi;
    Vector< std::unique_ptr< amrex::MultiFab > >  m_rhs;

    // Linear operator
    std::unique_ptr< amrex::MLNodeLaplacian > m_matrix;

    // Solver
    std::unique_ptr< amrex::MLMG > m_solver;


    void readParameters ();

    void computeRHS ( Vector< std::unique_ptr< amrex::MultiFab > >& a_vel,
                      Real a_time );

    void printInfo ();

};


#endif
