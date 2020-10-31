#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

class FieldRefinement : public RefinementCriteria::Register<FieldRefinement>
{
public:
    static std::string identifier() { return "FieldRefinement"; }

    FieldRefinement(const CFDSim& sim);

    virtual ~FieldRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_field;

    amrex::Vector<amrex::Real> m_field_error;
    amrex::Vector<amrex::Real> m_grad_error;

    int m_max_lev_field;
    int m_max_lev_grad;
};

} // namespace amr_wind

#endif /* FIELDREFINEMENT_H */
