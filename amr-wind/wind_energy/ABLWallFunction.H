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

    //! Return the plane-averaged computed friction velocity at any given
    //! instance
    amrex::Real utau() const { return m_utau; }

    //! Return the mean velocity used for friction velocity computations
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& umean() const
    {
        return m_umean;
    }

    //! Initialize the log-law height based on user inputs
    void init_log_law_height();

    //! Update the mean velocity at a given timestep
    void update_umean();

    /**
     * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
     * model consistent with Monin-Obukhov similarity theory. Wind Energy,
     * 20(3), 479–489. https://doi.org/10.1002/we.2017
     *
     * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
     * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
     * Demetri P Lalas and Corrado F Ratto, January 1996,
     * https://doi.org/10.1142/2975.
     */
    //! Return Monin-Obukhov Similarity function psi_m
    amrex::Real mo_psi_m(amrex::Real zeta);

    //! Return Monin-Obukhov Similarity function psi_h
    amrex::Real mo_psi_h(amrex::Real zeta);

    const amrex::FArrayBox& instplanar() const { return m_store_xy_vel_temp; }

    amrex::Real mean_windspeed() const { return m_mean_windspeed; }

    amrex::Real surface_temp_flux() const { return m_surf_temp_flux; }

    amrex::Real surface_temp() const { return m_surf_temp; }

    bool use_fch() const { return m_use_fch; }

    amrex::Real obukhov_length() const { return m_obukhov_length; }

private:
    void computeplanar();

    void computeusingheatflux();

    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! Mean velocity
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_umean{{0.0, 0.0, 0.0}};

    amrex::Real m_log_law_height{0.0}; ///< log-law height
    amrex::Real m_utau;                ///< Friction velocity

    amrex::Real m_kappa{0.41}; ///< von Karman constant
    amrex::Real m_z0{0.1};     ///< Roughness height

    int m_direction{2};   ///< Direction normal to wall
    bool m_use_fch{true}; ///< Use first cell height?

    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    int m_z_sample_index;
    amrex::Array<amrex::Real, 2> m_coeff_interp{{0.0, 0.0}};

    int m_ncells_x, m_ncells_y;

    amrex::Real m_mean_windspeed, m_mean_pot_temp;

    // amrex::Vector<amrex::Box> m_bx_z_sample;
    // amrex::Vector<amrex::FArrayBox>  m_store_xy_vel;
    // amrex::Vector<int> m_z_sample_index;

    amrex::Box m_bx_z_sample;
    amrex::FArrayBox m_store_xy_vel_temp;

    bool m_tempflux;
    amrex::Real m_surf_temp_flux;
    amrex::Real m_surf_temp_rate;
    amrex::Real m_surf_temp_rate_tstart;
    amrex::Real m_surf_temp_init;
    amrex::Real m_obukhov_length;
    amrex::Real m_psi_m, m_psi_h;
    amrex::Real m_surf_temp;
    amrex::Real m_ref_temp;

    amrex::Real m_gamma_m{5.0};
    amrex::Real m_gamma_h{5.0};
    amrex::Real m_beta_m{16.0};

    const int m_max_iter{25};
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
