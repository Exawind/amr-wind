#ifndef ISOPROBESAMPLER_H
#define ISOPROBESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind {
namespace sampling {

/** Sample data on a cloud of points
 *  \ingroup sampling
 *
 *  Reads a cloud of point locations from a text file and samples data on those
 *  points.
 */
class IsoProbeSampler : public SamplerBase::Register<IsoProbeSampler>
{
public:
    static const std::string identifier() { return "IsoProbeSampler"; }

    explicit IsoProbeSampler(const CFDSim&);

    virtual ~IsoProbeSampler();

    void initialize(const std::string& key) override;

    void sampling_locations(SampleLocType&) const override;
    void sampling_orientations(SampleLocType&) const override;

    void define_netcdf_metadata(const ncutils::NCGroup&) const override;

    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    int id() const override { return m_id; }
    int& id() override { return m_id; }
    int num_points() const override { return m_npts; }

private:
    SampleLocType m_probes, m_oris;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* ISOPROBESAMPLER_H */
