#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {

/** Wall-function computations for LES simulations of ABL
 *  \ingroup we_abl
 *
 *  This class performs the necessary computations at the beginning of
 *  predictor/corrector steps. The actual BC population in ghost cells is
 *  performed by ABLVelWallFunc BC interface class.
 */
class ABLWallFunction
{
public:
    ABLWallFunction(const CFDSim& sim);

    ~ABLWallFunction() = default;

    //! Return the height used to perform friction velocity computations
    amrex::Real log_law_height() const { return m_log_law_height; }

    //! Return the plane-averaged computed friction velocity at any given instance
    amrex::Real utau() const { return m_utau; }

    //! Return the mean velocity used for friction velocity computations
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& umean() const
    { return m_umean; }

    //! Initialize the log-law height based on user inputs
    void init_log_law_height();

    //! Update the mean velocity at a given timestep
    void update_umean(const FieldPlaneAveraging& pa);

private:
    const amrex::AmrCore& m_mesh;

    //! Mean velocity
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_umean{{0.0, 0.0, 0.0}};

    amrex::Real m_log_law_height{0.0}; ///< log-law height
    amrex::Real m_utau; ///< Friction velocity

    amrex::Real m_kappa{0.41}; ///< von Karman constant
    amrex::Real m_z0{0.1}; ///< Roughness height

    int m_direction{2}; ///< Direction normal to wall
    bool m_use_fch{true}; ///< Use first cell height?
};

/** Applies a shear-stress value at the domain boundary
 *  \ingroup field_bc we_abl
 *
 *  \sa ABLWallFunction
 */
class ABLVelWallFunc : public FieldBCIface
{
public:
    ABLVelWallFunc(Field& velocity, const ABLWallFunction& wall_func);

    void operator()(Field& velocity, const FieldState rho_state) override;

private:
    const ABLWallFunction& m_wall_func;
};

}

#endif /* ABLWALLFUNCTION_H */
