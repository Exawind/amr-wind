#ifndef DTUSPINNERSAMPLER_H
#define DTUSPINNERSAMPLER_H

#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::sampling {

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
    static std::string identifier() { return "DTUSpinnerSampler"; }

    explicit DTUSpinnerSampler(const CFDSim& /*sim*/);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    void sampling_locations(SampleLocType& /*locs*/) const override;

    /** Determine how the subsampling will be done
     *
     */
    void update_sampling_locations() override;

    int num_points() const override { return (m_npts * m_ns); }

    bool output_netcdf_field(
        double* /*unused*/, ncutils::NCVar& /*unused*/) override;

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/) const override;

private:
    int m_ns{1};
    amrex::Real m_dt_s{0.01};
    amrex::Real m_time_sampling{0};
};

} // namespace amr_wind::sampling

#endif /* DTUSPINNERSAMPLER_H */
