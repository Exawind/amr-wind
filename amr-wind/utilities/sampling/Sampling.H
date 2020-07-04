#ifndef SAMPLING_H
#define SAMPLING_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"

namespace amr_wind {
namespace sampling {

class Sampling : public PostProcessBase::Register<Sampling>
{
public:
    static const std::string identifier() { return "Sampling"; }

    Sampling(CFDSim&, const std::string&);

    virtual ~Sampling();

    void initialize();

    void post_advance_work();

protected:
    virtual void process_output();
    virtual void write_native();
    virtual void write_ascii();

private:
    CFDSim& m_sim;

    std::unique_ptr<SamplingContainer> m_scontainer;
    amrex::Vector<std::unique_ptr<SamplerBase>> m_samplers;

    amrex::Vector<std::string> m_var_names;
    amrex::Vector<Field*> m_fields;

    const std::string m_label;
    std::string m_out_fmt{"native"};

    int m_out_freq{100};
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLING_H */
