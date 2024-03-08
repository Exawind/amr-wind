#ifndef PROBESAMPLER_H
#define PROBESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind::sampling {

/** Sample data on a cloud of points
 *  \ingroup sampling
 *
 *  Reads a cloud of point locations from a text file and samples data on those
 *  points.
 */
class ProbeSampler : public SamplerBase::Register<ProbeSampler>
{
public:
    static std::string identifier() { return "ProbeSampler"; }

    //! Type of this sampling object
    std::string sampletype() const override { return identifier(); }

    explicit ProbeSampler(const CFDSim& /*unused*/);

    ~ProbeSampler() override;

    void initialize(const std::string& key) override;

    void sampling_locations(SampleLocType& /*locs*/) const override;

    void output_locations(SampleLocType& locs) const override
    {
        return sampling_locations(locs);
    }

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;

    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    long num_points() const override { return m_npts; }

    //! Number of output probe locations (after data reduction etc.)
    long num_output_points() const override { return m_npts; }

private:
    SampleLocType m_probes;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace amr_wind::sampling

#endif /* PROBESAMPLER_H */
