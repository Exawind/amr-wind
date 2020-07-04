#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <memory>

#include "amr-wind/core/Factory.H"

namespace amr_wind {

class CFDSim;

class PostProcessBase : public Factory<PostProcessBase, CFDSim&, const std::string&>
{
public:
    static const std::string base_identifier() { return "PostProcessBase"; }

    virtual ~PostProcessBase() = default;

    virtual void initialize() = 0;

    virtual void post_advance_work() = 0;
};

class PostProcessManager
{
public:
    PostProcessManager(CFDSim& sim);

    ~PostProcessManager() = default;

    void initialize();

    void post_advance_work();

private:
    CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<PostProcessBase>> m_post;
};

}

#endif /* POSTPROCESSING_H */
