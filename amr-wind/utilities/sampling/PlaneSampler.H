#ifndef PLANESAMPLER_H
#define PLANESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind {
namespace sampling {

class PlaneSampler : public SamplerBase::Register<PlaneSampler>
{
public:
    static const std::string identifier() { return "PlaneSampler"; }

    PlaneSampler(const CFDSim&);

    virtual ~PlaneSampler();

    void initialize(const std::string& key) override;

    void sampling_locations(SampleLocType&) const override;

    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    int id() const override { return m_id; }
    int& id() override { return m_id; }
    int num_points() const override { return m_npts; }

private:
    amrex::Vector<amrex::Real> m_axis1;
    amrex::Vector<amrex::Real> m_axis2;
    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<amrex::Real> m_normal{0.0, 0.0, 0.0};
    amrex::Vector<amrex::Real> m_poffsets;
    amrex::Vector<int> m_npts_dir;

    std::string m_label;

    int m_id{-1};
    int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* PLANESAMPLER_H */
