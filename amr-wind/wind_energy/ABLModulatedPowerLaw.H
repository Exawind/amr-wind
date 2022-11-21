#ifndef ABLMODULATEDPOWERLAW_H
#define ABLMODULATEDPOWERLAW_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_Gpu.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {

/** Interface for ABL modulated power law
 *  \ingroup we_abl
 *
 *  This class sets velocity and temperature at mass inflow bc's
 *
 */
class ABLModulatedPowerLaw
{

public:
    explicit ABLModulatedPowerLaw(CFDSim& /*sim*/);

    //! Execute initialization actions after mesh has been fully generated
    void post_init_actions();

    void pre_advance_work();

    void post_advance_work();

    void set_velocity(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

    void set_temperature(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

private:
    const CFDSim& m_sim;
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_temperature;

    amrex::Real m_zoffset{0.0};
    amrex::Real m_zref{90};
    amrex::Real m_shear_exp{0.1};
    //    amrex::Real m_umin{0.0};
    amrex::Real m_umax_factor{1.2};
    amrex::Real m_bulk_velocity{15.0};
    amrex::Real m_shearlayer_height{600.0};
    amrex::Real m_shearlayer_smear_thickness{30.0};
    amrex::Real m_wind_speed{8.0};
    amrex::Real m_wind_direction{270.0};
    amrex::Vector<amrex::Real> m_uvec{{8.0, 0.0, 0.0}};

    amrex::Real m_start_time{0.0};
    amrex::Real m_stop_time{std::numeric_limits<amrex::Real>::max()};
    amrex::Real m_degrees_per_sec{0.02};

    ///@{
    //! Lookup table containing temperature variation w.r.t. height above
    //! terrain
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    ///@}

    // Device copies of the above arrays
    amrex::Gpu::DeviceVector<amrex::Real> m_thht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_thvv_d;

    amrex::Real m_deltaT{0.8};
    amrex::Real m_theta_cutoff_height{250.0};
    amrex::Real m_theta_gauss_mean{0.0};
    amrex::Real m_theta_gauss_var{1.0};

    bool m_activate_mpl{false};
};

} // namespace amr_wind

#endif /* ABLMODULATEDPOWERLAW_H */
