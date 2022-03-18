#ifndef ABLFIELDINIT_H
#define ABLFIELDINIT_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize ABL fields
 */
class ABLFieldInit
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    ABLFieldInit(CFDSim&);

    void operator()(
        const int level,
        const amrex::MFIter& mfi,
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& density,
        const amrex::Array4<amrex::Real>& temperature) const;

    /** Add temperature perturbations
     *
     *  This uses amrex::Random and, therefore, executes the loop on CPU and
     *  pushes the field to device.
     */
    void perturb_temperature(
        const int lev, const amrex::Geometry& geom, Field& temperature) const;

    //! Flag indicating whether temperature field needs perturbations
    bool add_temperature_perturbations() const { return m_perturb_theta; }

    //! Initialize TKE field
    void init_tke(
        const int level,
        const amrex::Geometry& geom, 
	amrex::MultiFab& tke) const;

    //! Initialize SDR field
    void init_sdr(
        const int level,
        const amrex::Geometry& geom, 
        amrex::MultiFab& sdr) const;

    //! Initialize epsilon field
    void init_eps(
        const int level,
        const amrex::Geometry& geom, 
        amrex::MultiFab& sdr) const;

private:
    const FieldRepo& m_repo;

    //! Initial velocity components
    amrex::Vector<amrex::Real> m_vel;

    ///@{
    //! Lookup table containing temperature variation w.r.t. height above
    //! terrain
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    ///@}

    // Device copies of the above arrays
    amrex::Gpu::DeviceVector<amrex::Real> m_thht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_thvv_d;

    //! Using mesh mapping field or not
    bool m_mesh_mapping{false};

    //! Initial density field
    amrex::Real m_rho;

    //! Amplitude of fluctuations in x-direction
    amrex::Real m_deltaU{1.0};

    //! Amplitude of fluctuations in y-direction
    amrex::Real m_deltaV{1.0};

    //! Number of sinusoidal waves in x-direction
    amrex::Real m_Uperiods{4.0};

    //! Number of sinusoidal waves in y-direction
    amrex::Real m_Vperiods{4.0};

    //! Reference height for velocity perturbations
    amrex::Real m_ref_height{50.0};

    //! Amplitude of temperature perturbations
    amrex::Real m_deltaT{0.8};

    //! Mean for Gaussian number generator
    amrex::Real m_theta_gauss_mean{0.0};

    //! Variance for Gaussian number generator
    amrex::Real m_theta_gauss_var{1.0};

    //! Cutoff height for temperature fluctuations
    amrex::Real m_theta_cutoff_height{1.0e16};

    //! ABL initialization type
    std::string m_init_type{"const"};

    //! Initial constant value for kappa
    amrex::Real m_kappa{0.41};

    //! Initial constant value for kappa
    amrex::Real m_rough_z0{0.005};

    //! Initial forcing height
    amrex::Real m_force_height{32.1};

    //! Initial constant value for tke field
    amrex::Real m_tke_init{0.1};

    //! Initial constant value for eps field
    amrex::Real m_eps_init{0.1};

    //! Initial constant value for sdr field
    //! Default value set based on https://turbmodels.larc.nasa.gov/sst.html
    amrex::Real m_sdr_init{25.0};

    //! Initial value for ustar (used with geostrophic forcing only)
    amrex::Real m_init_prof_ustar{0.41};

    //! Initial value for C1 (see Parente2011)
    amrex::Real m_init_prof_c1{-0.04};

    //! Initial value for C1 (see Parente2011)
    amrex::Real m_init_prof_c2{0.52};

    //! Perturb initial velocity field with sinusoidal fluctuations
    bool m_perturb_vel{true};

    //! Perturb temperature field with random fluctuations
    bool m_perturb_theta{false};
};

} // namespace amr_wind

#endif /* ABLFIELDINIT_H */
