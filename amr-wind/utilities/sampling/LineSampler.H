#ifndef LINESAMPLER_H
#define LINESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind {
namespace sampling {

class LineSampler : public SamplerBase::Register<LineSampler>
{
public:
    static const std::string identifier() { return "LineSampler"; }

    LineSampler(const CFDSim&);

    virtual ~LineSampler();

    void initialize(const std::string& key) override;

    void sampling_locations(SampleLocType&) const override;

    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    int id() const override { return m_id; }
    int& id() override { return m_id; }
    int num_points() const override { return m_npts; }

private:
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

}
}

#endif /* LINESAMPLER_H */
