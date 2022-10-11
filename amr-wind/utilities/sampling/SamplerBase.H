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

    //! Name used to refer to this sampler (e.g., file I/O)
    [[nodiscard]] virtual std::string label() const = 0;
    virtual std::string& label() = 0;

    //! Unique numeric identifier used to track particles belonging to this
    //! instance
    [[nodiscard]] virtual int id() const = 0;
    virtual int& id() = 0;

    //! Total number of probes that belong to this sampling instance
    [[nodiscard]] virtual int num_points() const = 0;

    //! Read inputs and perform initialization actions
    virtual void initialize(const std::string& key) = 0;

    //! Populate the vector with coordinates of the sampling locations
    virtual void sampling_locations(SampleLocType&) const = 0;

    //! Update the sampling locations
    virtual void update_sampling_locations() {}

    //! Run specific output for the sampler
    virtual bool
    output_netcdf_field(double* /*unused*/, ncutils::NCVar& /*unused*/)
    {
        return true;
    }

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
