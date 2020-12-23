#ifndef ACTUATOR_TYPES_H
#define ACTUATOR_TYPES_H

#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_Gpu.H"
#include "AMReX_RealBox.H"

#include <set>

namespace amr_wind {

class CFDSim;

namespace actuator {

/** Abstract representation of an actuator type
 */
struct ActuatorType
{};

/** Abstract representation of a source type
 */
struct ActSrcType
{};

/** Actuator line representation
 */
struct ActSrcLine : ActSrcType
{
    static std::string identifier() { return "Line"; }

    //! Flag indicating if this is a line type
    static constexpr bool is_line = true;

    //! Flag indicating if this is a disk type
    static constexpr bool is_disk = false;
};

/** Acutator disk representation
 */
struct ActSrcDisk : ActSrcType
{
    static std::string identifier() { return "Disk"; }

    static constexpr bool is_line = true;
    static constexpr bool is_disk = false;
};

using VecList = amrex::Vector<amr_wind::vs::Vector>;
using TensorList = amrex::Vector<amr_wind::vs::Tensor>;
using DeviceVecList = amrex::Gpu::DeviceVector<amr_wind::vs::Vector>;
using DeviceTensorList = amrex::Gpu::DeviceVector<amr_wind::vs::Tensor>;

/** Actuator data at each node for an actuator component
 */
struct ActGrid
{
    //! Position vectors of the actuator forcing points on the grid
    VecList pos;

    //! Force vector at the forcing points
    VecList force;

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

    /** Convenience function to resize both force/velocity data to same locations
     */
    inline void resize(int num_pts) { resize(num_pts, num_pts); }
};

/** Basic information that is common to all actuator types
 */
struct ActInfo
{
    //! Unique string identifier for the actuator body (wing, turbine)
    std::string label;

    //! Set of MPI ranks where this actuator source is active
    std::set<int> procs;

    //! Bounding box used to determine processes where this turbine has
    //! influence
    amrex::RealBox bound_box;

    //! Unique integer identifier for the turbine
    int id{-1};

    //! Root process where this turbine is active
    int root_proc{-1};

    //! Flag indicating whether this actuator component is active in the current MPI rank
    bool actuator_in_proc{false};

    ActInfo(const std::string& label_in, const int id_in)
        : label(label_in), id(id_in)
    {}
};

/** Abstract representation of data holder for specific actuator types
 *
 */
template <typename ActTrait>
struct ActDataHolder
{
    //! Instance of the CFD simulation environment
    CFDSim& m_sim;

    //! Basic information about this actuator component in relation to the
    //! simulation (usually ActInfo)
    typename ActTrait::InfoType m_info;

    //! Nodal data on the actuator grid defined for this component (usually ActGrid)
    typename ActTrait::GridType m_grid;

    //! Additional data necessary for a given actuator type
    typename ActTrait::MetaType m_meta;

    /** Initialize the data structures
     *
     *  \param sim Instance of the CFD simulation environment
     *  \param label Unique string identifier for this actuator component
     *  \param id Unique global integer identifier for this actuator component
     */
    ActDataHolder(CFDSim& sim, const std::string& label, const int id)
        : m_sim(sim), m_info(label, id), m_grid(), m_meta()
    {}
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_TYPES_H */
