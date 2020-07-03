#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

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

}

#endif /* POSTPROCESSING_H */
