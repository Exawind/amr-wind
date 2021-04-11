#ifndef IBCONTAINER_H
#define IBCONTAINER_H

#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_AmrParticles.H"

namespace amr_wind {

class Field;

namespace ib {

class IB;

/** Cloud of IB geometry points.
 *
 *  \ingroup ibm
 *
 *  Holds the position vector and velocity information for a cloud of immersed
 *  boundary points belonging to all the immersed bodies within a given MPI
 * rank.
 */
struct IBCloud
{
    //! Number of nodes per turbine on this MPI rank
    amrex::Vector<int> num_pts;

    //! Global IDs of the ibms on this MPI rank
    amrex::Vector<int> global_id;

    //! Position vectors of all ibm nodes belonging to this MPI rank
    amrex::Vector<vs::Vector> position;

    //! Velocity vectors of all ibm nodes belonging to this MPI rank
    amrex::Vector<vs::Vector> velocity;

    //! Total number of ibms located on this MPI rank
    int num_objects;

    explicit ActuatorCloud(const int nobjects);
};

//! Number or real entries in Array of Structs (AOS)
static constexpr int NumPStructReal = AMREX_SPACEDIM;
//! Number of integer entries in Array of Structs (AOS)
static constexpr int NumPStructInt = 1;
//! Number of real entries in Struct of Arrays (SOA)
static constexpr int NumPArrayReal = 0;
//! Number of int entries in Struct of Arrays (SOA)
static constexpr int NumPArrayInt = 0;

/** Specialization of AmrParticleContainer for sampling velocities.
 *
 *  \ingroup actuator
 */
class IBContainer
    : public amrex::AmrParticleContainer<
          NumPStructReal,
          NumPStructInt,
          NumPArrayReal,
          NumPArrayInt>
{
public:
    friend class IB;

    explicit IBContainer(amrex::AmrCore& mesh, const int num_objects);

    void post_regrid_actions();

    void initialize_container();

    void reset_container();

    void update_positions();

    void sample_velocities(const Field& vel);

    int num_actuator_points() const { return m_data.position.size(); }

    // public for CUDA, not safe for general access

    void interpolate_velocities(const Field& vel);

    void populate_vel_buffer();

    void initialize_particles(const int total_pts);

protected:
    void compute_local_coordinates();

    // Accessor to allow unit testing
    IBCloud& point_data() { return m_data; }

private:
    amrex::AmrCore& m_mesh;

    // Object that holds the position and velocity information
    ActuatorCloud m_data;

    //! Coordinates of points that are contained within boxes assigned to each
    //! MPI rank
    amrex::Vector<vs::Vector> m_proc_pos;

    //! Device view of the process position vectors
    amrex::Gpu::DeviceVector<vs::Vector> m_pos_device;

    amrex::Vector<int> m_proc_offsets;
    amrex::Gpu::DeviceVector<int> m_proc_offsets_device;

    //! Flag indicating whether memory has allocated for all data structures
    bool m_container_initialized{false};

    //! Flag indicating whether the particles are scattered throughout the
    //! domain, or if they have been recalled to the original MPI rank
    bool m_is_scattered{false};
};

} // namespace ib
} // namespace amr_wind

#endif /* ACTUATORCONTAINER_H */
