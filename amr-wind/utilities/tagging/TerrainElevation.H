/*
 * NOTE: This class is a work in progress.
 * In the future, TerrainElevation may be extended to handle additional
 * terrain-related data such as roughness, and could be used by other components
 * of the codebase, including TerrainDrag, MetMastForcing, and similar modules.
 */

#ifndef TERRAINELEVATION_H_
#define TERRAINELEVATION_H_

#include "AMReX.H"
#include <AMReX_Vector.H>
#include <memory>
#include <string>

namespace amr_wind {

class TerrainElevation
{
public:
    TerrainElevation() = default;

    void load_from_file(const std::string& filename);
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool is_loaded() const
    {
        return m_loaded;
    }

    static std::shared_ptr<TerrainElevation>& get_instance();
    static void ensure_loaded(const std::string& filename);

    // Accessors (const and non-const)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const amrex::Vector<amrex::Real>&
    x() const
    {
        return m_x;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const amrex::Vector<amrex::Real>&
    y() const
    {
        return m_y;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const amrex::Vector<amrex::Real>&
    z() const
    {
        return m_z;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Vector<amrex::Real>& x()
    {
        return m_x;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Vector<amrex::Real>& y()
    {
        return m_y;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Vector<amrex::Real>& z()
    {
        return m_z;
    }

private:
    amrex::Vector<amrex::Real> m_x;
    amrex::Vector<amrex::Real> m_y;
    amrex::Vector<amrex::Real> m_z;
    bool m_loaded{false};
};

// Singleton accessor
std::shared_ptr<TerrainElevation>& get_terrain_elevation();

} // namespace amr_wind

#endif /* TERRAINELEVATION_H_ */