#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/MOData.h"

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

    MOData& mo() { return m_mo; }
    const MOData& mo() const { return m_mo; }

    //! Return the plane-averaged computed friction velocity at any given
    //! instance
    amrex::Real utau() const { return m_mo.utau; }

    //! Initialize the log-law height based on user inputs
    void init_log_law_height();

    //! Update the mean velocity at a given timestep
    void update_umean(const VelPlaneAveraging& vpa, const FieldPlaneAveraging& tpa);

    void computeplanar();

    void computeusingheatflux();

    const amrex::FArrayBox& instplanar() const { return m_store_xy_vel_temp; }

private:

    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! Monin-Obukhov instance
    MOData m_mo;

    int m_direction{2};   ///< Direction normal to wall
    bool m_use_fch{true}; ///< Use first cell height?

    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    int m_z_sample_index;
    amrex::Array<amrex::Real, 2> m_coeff_interp{{0.0, 0.0}};

    int m_ncells_x, m_ncells_y;

    amrex::Box m_bx_z_sample;
    amrex::FArrayBox m_store_xy_vel_temp;

    bool m_tempflux;
    amrex::Real m_surf_temp_rate;
    amrex::Real m_surf_temp_rate_tstart;
    amrex::Real m_surf_temp_init;
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

class ABLTempWallFunc : public FieldBCIface
{
public:
    ABLTempWallFunc(Field& temperature, const ABLWallFunction& wall_func);

    void operator()(Field& temperature, const FieldState rho_state) override;

private:
    const ABLWallFunction& m_wall_func;
};

} // namespace amr_wind

#endif /* ABLWALLFUNCTION_H */
