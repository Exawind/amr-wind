#ifndef VOLUMESAMPLER_H
#define VOLUMESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind::sampling {

/** Sample data on uniform grids
 *  \ingroup sampling
 *
 *  Defines probe locations on a volume defined by an
 *  ` volume that originates from the point `lo` and extends to `hi`.
 *  The resolution of the volume is defined by the vector npts_dir.
 *
 */
class VolumeSampler : public SamplerBase::Register<VolumeSampler>
{
public:
    static std::string identifier() { return "VolumeSampler"; }

    explicit VolumeSampler(const CFDSim& /*unused*/);

    ~VolumeSampler() override;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*locs*/) const override;

    void output_locations(SampleLocType& locs) const override
    {
        return sampling_locations(locs);
    }

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;

    //! Name of this sampling object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    std::string sampletype() const override { return identifier(); }

    //! Unique identifier for this set of probe locations
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of probe locations along the line
    long num_points() const override { return m_npts; }

    //! Number of output probe locations (after data reduction etc.)
    long num_output_points() const override { return m_npts; }

private:
    amrex::Vector<amrex::Real> m_hi;
    amrex::Vector<amrex::Real> m_lo;
    amrex::Vector<int> m_npts_dir;

    std::string m_label;

    int m_id{-1};
    int m_npts{0};
};

} // namespace amr_wind::sampling

#endif /* VOLUMESAMPLER_H */
