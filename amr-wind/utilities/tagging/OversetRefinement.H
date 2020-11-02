#ifndef OVERSETREFINEMENT_H
#define OVERSETREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {

class OversetRefinement : public RefinementCriteria::Register<OversetRefinement>
{
public:
    static std::string identifier() { return "OversetRefinement"; }

    OversetRefinement(const CFDSim& sim);

    virtual ~OversetRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    int m_max_lev;
};

} // namespace amr_wind

#endif /* OVERSETREFINEMENT_H */
