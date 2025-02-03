#ifndef OCEANWAVESBOUNDARY_H
#define OCEANWAVESBOUNDARY_H

#include "amr-wind/core/Field.H"
#include "amr-wind/core/IntField.H"
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

    void post_init_actions();

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

    void set_inflow_sibling_velocity(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs) const;

    void record_boundary_data_time(const amrex::Real time)
    {
        m_bndry_time = time;
    }

    amrex::Real get_boundary_data_time() const { return m_bndry_time; }

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_ow_velocity;
    Field& m_ow_vof;
    IntField* m_terrain_blank_ptr{nullptr};

    // Should be active unless boundary planes are used
    bool m_activate_ow_bndry{true};

    // Check for single-phase simulation
    bool m_vof_exists{true};

    // Check for terrain
    bool m_terrain_exists{false};

    // Store time corresponding to current boundary data
    amrex::Real m_bndry_time{0.0};

    // Will only remain negative if not needed
    amrex::Real m_rho1{-1.0};
};

} // namespace amr_wind

#endif /* OCEANWAVESBOUNDARY_H */
