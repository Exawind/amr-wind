#ifndef TRACERREFINEMENT_H
#define TRACERREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {

class TracerRefinement : public RefinementCriteria::Register<TracerRefinement>
{
public:
    static std::string identifier() { return "TracerRefinement"; }

    explicit TracerRefinement(const CFDSim& sim);

    virtual ~TracerRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    int m_max_lev;
    amrex::Real m_value{0.1};
    Field& m_passive_scalar;
};

} // namespace amr_wind

#endif /* TRACERREFINEMENT_H */
