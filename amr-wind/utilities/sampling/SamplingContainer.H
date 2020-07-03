#ifndef SAMPLINGCONTAINER_H
#define SAMPLINGCONTAINER_H

#include <memory>

#include "AMReX_AmrParticles.H"

namespace amr_wind {

class Field;

namespace sampling {

class SamplerBase;

static constexpr int SNStructReal = 0;
static constexpr int SNStructInt = 2;
static constexpr int SNArrayReal = 0;
static constexpr int SNArrayInt = 0;

struct IIx
{
    enum Indices {
        sid = 0, ///< Identifier of the set this particle belongs to
        nid      ///< Index within the set for this particle
    };
};

class SamplingContainer
    : public amrex::AmrParticleContainer<
          SNStructReal,
          SNStructInt,
          SNArrayReal,
          SNArrayInt>
{
public:
    SamplingContainer(amrex::AmrCore& mesh)
        : AmrParticleContainer<
              SNStructReal,
              SNStructInt,
              SNArrayReal,
              SNArrayInt>(&mesh)
        , m_mesh(mesh)
    {}

    /** Define the container and allocate memory for runtime components
     */
    void setup_container(const int num_real_components,
                         const int num_int_components = 0);

    /** Create particle information for all the sampling locations
     */
    void initialize_particles(
        const amrex::Vector<std::unique_ptr<SamplerBase>>&);

    //! Perform field interpolation to sampling locations
    void interpolate_fields(const amrex::Vector<Field*> fields);

    int num_sampling_particles() const { return m_total_particles; }

private:
    amrex::AmrCore& m_mesh;

    int m_total_particles;
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLINGCONTAINER_H */
