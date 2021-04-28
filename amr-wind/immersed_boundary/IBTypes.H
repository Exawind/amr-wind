#ifndef IBTYPES_H
#define IBTYPES_H

#include "amr-wind/core/Slice.H"
#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_Gpu.H"
#include "AMReX_RealBox.H"

#include <set>

namespace amr_wind {

class CFDSim;

namespace ib {

/** Abstract representation of an immersed boundary type.
 *
 *  \ingroup immersed boundary
 */
struct IBType
{};

// FIXME: Does immersed boundary need any notion of multiple source types ?
/** Abstract representation of a source type.
 *
 *  \ingroup immersed boundary
 */
struct IBSrcType
{};

/** Immersed boundary representation.
 *
 *  \ingroup Immersed boundary
 */
struct IBSrcLine : IBSrcType
{
    // FIXME: Change identifier name
    static std::string identifier() { return "Line"; }
};

using RealList = amrex::Vector<amrex::Real>;
using RealSlice = ::amr_wind::utils::Slice<amrex::Real>;
using VecList = amrex::Vector<amr_wind::vs::Vector>;
using VecSlice = ::amr_wind::utils::Slice<amr_wind::vs::Vector>;
using TensorList = amrex::Vector<amr_wind::vs::Tensor>;
using TensorSlice = ::amr_wind::utils::Slice<amr_wind::vs::Tensor>;
using DeviceVecList = amrex::Gpu::DeviceVector<amr_wind::vs::Vector>;
using DeviceTensorList = amrex::Gpu::DeviceVector<amr_wind::vs::Tensor>;

/** Immersed boundary data at each node for an immersed boundary component.
 *
 *  \ingroup immersed boundary
 */
struct IBGrid
{
    //! Position vectors of the immersed boundary forcing points on the grid
    VecList pos;

    //! Force vector at the forcing points
    VecList force;

    // FIXME: Does immersed boundary need different such parameters ?
    //! Gaussian smearing parameter
    VecList epsilon;

    //! Transformation matrix at the forcing points
    TensorList orientation;

    //! Position vectors for the points where velocity is sampled
    VecList vel_pos;

    //! Velocity vector at the sampled locations
    VecList vel;

    /** Helper method to resize the data arrays defined on the grid
     *
     *  \params Number of nodes that contain forcing data
     *  \params Number of nodes where velocity field is sampled
     */
    void resize(int num_force_pts, int num_vel_pts)
    {
        pos.resize(num_force_pts);
        force.resize(num_force_pts);
        epsilon.resize(num_force_pts);
        orientation.resize(num_force_pts);
        vel_pos.resize(num_vel_pts);
        vel.resize(num_vel_pts);
    }

    /** Convenience function to resize both force/velocity data to same
     * locations
     */
    inline void resize(int num_pts) { resize(num_pts, num_pts); }
};

/** Basic information that is common to all immersed boundary types.
 *
 *  \ingroup immersed boundary
 */
struct IBInfo
{
    //! Unique string identifier for the immersed boundary body (wing, turbine)
    std::string label;

    //! Set of MPI ranks where this immersed boundary source is active
    std::set<int> procs;

    //! Bounding box used to determine processes where this turbine has
    //! influence
    amrex::RealBox bound_box;

    //! Unique integer identifier for the turbine
    int id{-1};

    //! Root process where this turbine is active
    int root_proc{-1};

    //! Flag indicating whether this is root proc
    bool is_root_proc{false};

    //! Flag indicating whether this immersed boundary component is active in the current
    //! MPI rank
    bool ib_in_proc{false};

    //! Flag indicating whether this process requires velocities sampled at
    //! immersed boundary point
    bool sample_vel_in_proc{false};

    IBInfo(const std::string& label_in, const int id_in)
        : label(label_in), id(id_in)
    {}
};

/** Abstract representation of data holder for specific immersed boundary types.
 *
 *  \ingroup immersed boundary
 *
 */
template <typename IBTrait>
class IBDataHolder
{
private:
    //! Instance of the CFD simulation environment
    CFDSim& m_sim;

    //! Basic information about this immersed boundary component in relation to the
    //! simulation (usually IBInfo)
    typename IBTrait::InfoType m_info;

    //! Nodal data on the immersed boundary grid defined for this component (usually
    //! ActGrid)
    typename IBTrait::GridType m_grid;

    //! Additional data necessary for a given actuator type
    typename IBTrait::MetaType m_meta;

public:
    /** Initialize the data structures
     *
     *  \param sim Instance of the CFD simulation environment
     *  \param label Unique string identifier for this immersed boundary component
     *  \param id Unique global integer identifier for this immersed boundary component
     */
    IBDataHolder(CFDSim& sim, const std::string& label, const int id)
        : m_sim(sim), m_info(label, id), m_grid(), m_meta()
    {}

    CFDSim& sim() { return m_sim; }
    const CFDSim& sim() const { return m_sim; }

    typename IBTrait::InfoType& info() { return m_info; }
    const typename IBTrait::InfoType& info() const { return m_info; }

    typename IBTrait::GridType& grid() { return m_grid; }
    const typename IBTrait::GridType& grid() const { return m_grid; }

    typename IBTrait::MetaType& meta() { return m_meta; }
    const typename IBTrait::MetaType& meta() const { return m_meta; }
};

} // namespace ib
} // namespace amr_wind

#endif /* IBTYPES_H */
