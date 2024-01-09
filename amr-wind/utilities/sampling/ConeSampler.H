#ifndef CONESAMPLER_H
#define CONESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segments as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class ConeSampler : public SamplerBase::Register<ConeSampler>
{
public:
    static std::string identifier() { return "ConeSampler"; }

    explicit ConeSampler(const CFDSim& /*sim*/);

    ~ConeSampler() override;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    void check_bounds();

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
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/) const override;

    //! Name of this sampling object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    std::string sampletype() const override { return identifier(); }

    //! Unique identifier for this set of probe locations
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of probe locations
    long num_points() const override
    {
        return m_npts + m_ntheta * (m_npts - 1) * (m_nphi - 1);
    }

    //! Number of probe locations on cone axis
    long num_points_axis() const { return m_npts; }

    //! Number of output probe locations (after data reduction etc.)
    long num_output_points() const override { return m_npts; }

protected:
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    std::string m_label;
    int m_id{-1};

    int m_npts{0};
    int m_ntheta{0};
    int m_nphi{0};
    double m_phi{0.0};
};

} // namespace amr_wind::sampling

#endif /* CONESAMPLER_H */
