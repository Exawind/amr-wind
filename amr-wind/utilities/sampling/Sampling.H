#ifndef SAMPLING_H
#define SAMPLING_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"

namespace amr_wind {
namespace sampling {

/** Collection of data sampling objects
 *
 *  A concrete impelementation of the post-processing interface that deals with
 *  data probes. Multiple instances of sampling entities (e.g., point cloud or
 *  probes, line, or plane samplers) can be registered. Internally, the data
 *  probles are represented as amr_wind::sampling::SamplingContainer for
 *  tracking their position within the CFD mesh and updating the fields of
 *  interest.
 */
class Sampling : public PostProcessBase::Register<Sampling>
{
public:
    static const std::string identifier() { return "Sampling"; }

    Sampling(CFDSim&, const std::string&);

    virtual ~Sampling();

    //! Read user inputs and create the different data probe instances
    void initialize();

    //! Interpolate fields at a given timestep and output to disk
    void post_advance_work();

    // Public for CUDA

    //! Write sampled data in binary format
    virtual void impl_write_native();

protected:
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

    //! Number of particles:
    size_t m_total_particles;
    //! Frequency of data sampling and output
    int m_out_freq{100};
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLING_H */
