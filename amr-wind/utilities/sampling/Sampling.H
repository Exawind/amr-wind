#ifndef SAMPLING_H
#define SAMPLING_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"

/**
 *  \defgroup sampling Data-sampling utilities
 *  Data-sampling utilities
 *
 *  Sampling contains various data-probes sampling utilities that are
 *  implemented in using the [AMReX
 *  particles](https://amrex-codes.github.io/amrex/docs_html/Particle.html)
 *  interface. It supports output in native AMReX particles binary format as
 *  well as NetCDF format.
 *
 *  \ingroup utilities
 */

namespace amr_wind::sampling {

/** Collection of data sampling objects
 *  \ingroup sampling
 *
 *  A concrete impelementation of the post-processing interface that deals with
 *  data probes. Multiple instances of sampling entities (e.g., point cloud or
 *  probes, line, or plane samplers) can be registered. Internally, the data
 *  probles are represented as amr_wind::sampling::SamplingContainer for
 *  tracking their position within the CFD mesh and updating the fields of
 *  interest.
 *
 *  This class manages the interpolation and the field management. The actual
 *  definition of probes are managed by subclasses of SamplerBase.
 */
class Sampling : public PostProcessBase::Register<Sampling>
{
public:
    static std::string identifier() { return "Sampling"; }

    Sampling(CFDSim& /*sim*/, std::string /*label*/);

    ~Sampling() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    //! Interpolate fields at a given timestep and output to disk
    void post_advance_work() override;

    //! Actions to perform post regrid e.g. redistribute particles
    void post_regrid_actions() override;

    // Public for CUDA

    //! Write sampled data in binary format
    virtual void impl_write_native();

    int num_total_particles() const { return m_total_particles; }

    const amrex::Vector<std::string>& var_names() const { return m_var_names; }

protected:
    //! Update the container by re-initializing the particles
    void update_container();

    //! Conduct work to update the particles
    void update_sampling_locations();

    //! Output data based on user-defined format
    virtual void process_output();

    //! Prepare NetCDF metadata
    virtual void prepare_netcdf_file();

    //! Write sampled data into a NetCDF file
    void write_netcdf();

    /** Output sampled data in ASCII format
     *
     *  Note that this should be used for debugging only and not in production
     *  runs as it can have significant impacts on code performance.
     */
    virtual void write_ascii();

    SamplingContainer& sampling_container() { return *m_scontainer; }

private:
    CFDSim& m_sim;

    std::unique_ptr<SamplingContainer> m_scontainer;
    amrex::Vector<std::unique_ptr<SamplerBase>> m_samplers;

    //! List of variable names for output
    amrex::Vector<std::string> m_var_names;

    //! List of fields to be sampled for this collection of probes
    amrex::Vector<Field*> m_fields;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! Format of the data output (native, ascii, netcdf, etc.)
#ifdef AMR_WIND_USE_NETCDF
    std::string m_out_fmt{"netcdf"};
    std::string m_ncfile_name;
#else
    std::string m_out_fmt{"native"};
#endif

    // number of field components
    int m_ncomp{0};

    //! Number of particles:
    size_t m_total_particles{0};

    //! Frequency of data sampling and output
    int m_out_freq{100};
};

} // namespace amr_wind::sampling

#endif /* SAMPLING_H */
