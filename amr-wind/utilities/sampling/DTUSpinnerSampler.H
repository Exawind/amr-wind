#ifndef DTUSPINNERSAMPLER_H
#define DTUSPINNERSAMPLER_H

//~ #include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace sampling {

/** Sample data along a line that moves in time
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by a time history from `m_azimuth_table`
 *  and `m_elevation_table` with a length of `m_length`. This determines the
 * `start` and `end` points at every time-step and divides the line into equal
 * segements as determined by `num_points`. The length of each segment is given
 * by \f[ m_length / (num_points - 1) \f]
 */
class DTUSpinnerSampler : public LidarSampler
{
public:
    static const std::string identifier() { return "DTUSpinnerSampler"; }

    explicit DTUSpinnerSampler(const CFDSim&);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    //~ void update_sampling_locations() override;
    void sampling_locations(SampleLocType&) const override;

    /** Determine how the subsampling will be done
     *
     */
    void subsampling() override;
    void update_sampling_locations() override;

    int num_points() const override { return (m_npts * m_ns); }

    void define_netcdf_metadata(const ncutils::NCGroup&) const override;
    void populate_netcdf_metadata(const ncutils::NCGroup&) const override;
    void
    output_netcdf_data(const ncutils::NCGroup&, const size_t) const override;

private:
    int m_ns{1};
    amrex::Real m_dt_s{0.01};
    amrex::Real m_time_sampling{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* DTUSPINNERSAMPLER_H */
