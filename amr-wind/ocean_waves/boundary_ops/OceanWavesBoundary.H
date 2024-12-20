#ifndef OCEANWAVESBOUNDARY_H
#define OCEANWAVESBOUNDARY_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_Gpu.H"
#include <AMReX_BndryRegister.H>

namespace amr_wind {

/** Interface for Ocean Waves boundary data
 *
 *  This class performs the necessary field operations to provide
 *  available wave data at the boundaries
 *
 */
class OceanWavesBoundary
{

public:
    explicit OceanWavesBoundary(CFDSim& /*sim*/);

    //! Execute initialization actions after mesh has been fully generated
    void post_init_actions();

    void pre_advance_work();

    void pre_predictor_work();

    void post_advance_work();

    void initialize_data();

    void set_velocity(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        const int dcomp = 0,
        const int orig_comp = 0) const;

    void set_vof(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

    void set_density(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_ow_velocity;
    Field& m_ow_vof;

    // Should be active unless boundary planes are used
    bool m_activate_ow_bndry{true};

    // Will only remain negative if not needed
    amrex::Real m_rho1{-1.0};
};

} // namespace amr_wind

#endif /* OCEANWAVESBOUNDARY_H */
