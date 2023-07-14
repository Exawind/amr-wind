#ifndef SAMPLINGCONTAINER_H
#define SAMPLINGCONTAINER_H

#include <memory>

#include "AMReX_AmrParticles.H"

namespace amr_wind {

class Field;

namespace sampling {

class SamplerBase;

static constexpr int SNStructReal = 0;
static constexpr int SNStructInt = 3;
static constexpr int SNArrayReal = 0;
static constexpr int SNArrayInt = 0;

/** Indexing for the integer variables defined on the AMReX particle
 *  \ingroup sampling
 */
struct IIx
{
    enum Indices {
        uid = 0, ///< Unique identifier for this particle
        sid,     ///< Identifier of the set this particle belongs to
        nid      ///< Index within the set for this particle
    };
};

/** Specialization of the AMReX ParticleContainer object for sampling data
 *  \ingroup sampling
 *
 *  Represents the data probes as a collection of particles that are uniquely
 *  identified by three integer identifiers:
 *
 *   - a universal identifier (UID) that is unique to each particle,
 *   - a set identifier that represents a group (e.g., line, plane) of particles
 *   - a node identifier that represents indexing of this particle within the
 * group
 *
 *  The quantities of interest are themselves stored as runtime components based
 *  on user-defined inputs.
 *
 *
 *  Notes:
 *
 *   - The implementation uses linear interpolation in three directions to
 *     determine the data at a given probe location
 *
 *   - For non-nodal fields, the current implementation requires at-least one
 *     ghost cell to allow linear interpolation.
 *
 *   - Interpolation near domain boundaries does not currently handle `hoextrap`
 */
class SamplingContainer
    : public amrex::AmrParticleContainer<
          SNStructReal,
          SNStructInt,
          SNArrayReal,
          SNArrayInt>
{
public:
    explicit SamplingContainer(amrex::AmrCore& mesh)
        : amrex::AmrParticleContainer<
              SNStructReal,
              SNStructInt,
              SNArrayReal,
              SNArrayInt>(&mesh)
        , m_mesh(mesh)
    {}

    /** Define the container and allocate memory for runtime components
     */
    void setup_container(
        const int num_real_components, const int num_int_components = 0);

    /** Create particle information for all the sampling locations
     */
    void initialize_particles(
        const amrex::Vector<std::unique_ptr<SamplerBase>>& /*samplers*/);

    //! Perform field interpolation to sampling locations
    void interpolate_fields(const amrex::Vector<Field*> fields);

    //! Populate the buffer with data for all the particles
    void populate_buffer(std::vector<double>& buf);

    int num_sampling_particles() const { return m_total_particles; }

    int& num_sampling_particles() { return m_total_particles; }

private:
    amrex::AmrCore& m_mesh;

    int m_total_particles{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLINGCONTAINER_H */
