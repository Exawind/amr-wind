#ifndef PDEBASE_H
#define PDEBASE_H

#include <string>

#include "Factory.H"
#include "incflo_enums.H"
#include "FieldDescTypes.H"
#include "FieldRepo.H"
#include "PDEHelpers.H"
#include "SimTime.H"
#include "CollMgr.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Abstract interface defining a PDE and its solution process withn AMR-Wind
 *
 *  The exact behavior is implemented by subclass PDESystem that is specialized
 *  using two traits: the PDE type trait, and the numerical scheme trait.
 *  Currently, two numerical scheme traits are supported: Godunov and MOL.
 *
 *  Each PDESystem instance contains its own PDEFields object that holds fields
 *  unique to the PDE. These fields are: the variable solved by the PDE (e.g.,
 *  velocity, temperature, etc.), the effective dynamic viscosity, the source
 *  terms, the diffusion term, and the convective term. A PDE may declare
 *  additional fields that are not stored in this object, but can be queried
 *  using the FieldRepo instance. For example, the incompressible Navier-Stokes
 *  (ICNS) always declars density and pressure in addition to the velocity
 *  fields.
 *
 */
class PDEBase : public Factory<PDEBase, CFDSim&, const int>
{
public:
    virtual ~PDEBase() = default;

    //! Return the object that holds references to fields for this PDE system
    virtual PDEFields& fields() = 0;

    //! Perform initialization actions for a PDE after the mesh is generated
    virtual void initialize() = 0;

    //! Preform updates specific to PDE after a regrid is performed
    virtual void post_regrid_actions() = 0;

    //! Compute the source (forcing) terms for this PDE
    virtual void compute_source_term(const FieldState fstate) = 0;

    //! Compute the effective dynamic viscosity for this PDE
    virtual void compute_mueff(const FieldState fstate) = 0;

    //! Compute the diffusion term used in the RHS of the PDE system
    virtual void compute_diffusion_term(const FieldState fstate) = 0;

    //! Compute the time derivative and advective term for the PDE system
    virtual void compute_advection_term(const FieldState fstate) = 0;

    /** Combine the source, diffusion, and advection term to obtain RHS
     *
     *  This method behaves differently depending upon whether the diffusion
     *  term is treated implicitly, with Crank-Nicolson scheme, or explicitly.
     *  The predictor step only uses the "old" states for RHS computation.
     */
    virtual void compute_predictor_rhs(const DiffusionType difftype) = 0;

    /** Combine the source, diffusion, and advection term to obtain RHS
     *
     *  This method behaves differently depending upon whether the diffusion
     *  term is treated implicitly, with Crank-Nicolson scheme, or explicitly.
     *  The corrector step uses both "old" state as well as the updated "new"
     *  state from predictor step to compute the RHS term.
     */
    virtual void compute_corrector_rhs(const DiffusionType difftype) = 0;

    //! Solve the diffusion linear system and update the field
    virtual void solve(const amrex::Real dt) = 0;

    //! Base class identifier used for factory registration interface
    static std::string base_identifier() { return "PDESystem"; }
};

/** PDE systems manager
 *
 *  This class manages a collection of PDEs and provides an interface to
 *  register PDEs. An incompressible N-S equation is always registered and
 *  treated as a special case. The remaining PDEs are created on-demand and are
 *  accessed through the equation systems vector.
 */
class PDEMgr : public CollMgr<PDEMgr, PDEBase>
{
public:
    PDEMgr(CFDSim& sim);

    ~PDEMgr() = default;

    //! Return the incompressible Navier-Stokes instance
    PDEBase& icns() { return *m_icns; }

    //! Create the incompressible Navier-Stokes instance and return its reference
    PDEBase& register_icns();

    //! Register a new PDE instance and return its reference
    PDEBase& register_transport_pde(const std::string& pde_name);

    /** Return the vector containing all registered PDE instances
     *
     *  Note that this does not contain the ICNS system as it is treated
     *  separately for algorithmic reasons.
     */
    TypeVector& scalar_eqns() { return m_obj_vec; }

private:
    //! Instance of the CFD simulation controller
    CFDSim& m_sim;

    //! ICNS instance
    std::unique_ptr<PDEBase> m_icns;

    //! Unique identifier for the numerical scheme
    std::string m_scheme;

    //! Unique integer for the problem type (to be deprecated)
    int m_probtype;

    //! Flag indicating whether Godunov scheme is active
    bool m_use_godunov{false};
};

} // namespace pde
} // namespace amr_wind

#endif /* PDEBASE_H */
