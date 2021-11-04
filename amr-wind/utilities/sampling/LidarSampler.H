#ifndef LIDARSAMPLER_H
#define LIDARSAMPLER_H

#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segements as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class LidarSampler : public LineSampler
{
public:
    static const std::string identifier() { return "LidarSampler"; }

    explicit LidarSampler(const CFDSim&);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    void update_sampling_locations() override;

private:
    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<amrex::Real> m_time_table;
    amrex::Vector<amrex::Real> m_azimuth_table;
    amrex::Vector<amrex::Real> m_elevation_table;
    amrex::Real m_length{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* LIDARSAMPLER_H */
