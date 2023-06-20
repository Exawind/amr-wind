#ifndef SAMPLERBASE_H
#define SAMPLERBASE_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

namespace amr_wind {

class CFDSim;

namespace sampling {

/** Abstract representation of data probes to sample flow data
 *  \ingroup sampling
 *
 *  This class defines the interface by which the positions of the data-probes
 *  are defined and access from the Sampling class. All data probe definitions
 *  must be sub-classes of this base class.
 */
class SamplerBase : public Factory<SamplerBase, CFDSim&>
{
public:
    using SampleLocType =
        amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>>;

    static std::string base_identifier() { return "SamplerBase"; }

    ~SamplerBase() override = default;

    //! Name used to refer to this sampler (e.g., myline1)
    virtual std::string label() const = 0;
    virtual std::string& label() = 0;

    //! Class name of this sampler (e.g., LineSampler)
    virtual std::string sampletype() const = 0;

    //! Unique numeric identifier used to track particles belonging to this
    //! instance
    virtual int id() const = 0;
    virtual int& id() = 0;

    //! Total number of probes that belong to this sampling instance
    virtual int num_points() const = 0;

    //! Total number of probes to be output
    //! Note: Same as num_points() for sampling without data reduction
    virtual int num_output_points() const = 0;

    //! Read inputs and perform initialization actions
    virtual void initialize(const std::string& key) = 0;

    //! Populate the vector with coordinates of the sampling locations
    virtual void sampling_locations(SampleLocType&) const = 0;

    //! Populate the vector with coordinates of the output locations
    virtual void output_locations(SampleLocType&) const = 0;

    //! Update the sampling locations
    virtual void update_sampling_locations() {}

    //! Run specific output for the sampler
    virtual bool output_netcdf_field(
        const std::vector<double>& /*unused*/,
        ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/)
    {
        return true;
    }

    //! Run data modification for specific sampler?
    virtual bool do_data_modification() { return false; }
    virtual bool do_convert_velocity_los() { return false; }

    //! Sample buffer modification instructions
    virtual std::vector<double> modify_sample_data(
        const std::vector<double>& sampledata, const std::string& /*unused*/)
    {
        return sampledata;
    }

    virtual void
    calc_lineofsight_velocity(const std::vector<std::vector<double>>&)
    {}

    //! Populate metadata in the NetCDF file
    virtual void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const
    {}
    virtual void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const
    {}
    virtual void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
    {}
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLERBASE_H */
