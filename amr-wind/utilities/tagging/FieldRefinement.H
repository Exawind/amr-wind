#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = FieldRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field value
 *  # tagging/t1.field_error = 10.0 10.0

 *  # Tagging based on local gradient
 *  tagging/t1.grad_error = 0.1 0.1 0.1
 *  ```
 */
class FieldRefinement : public RefinementCriteria::Register<FieldRefinement>
{
public:
    static std::string identifier() { return "FieldRefinement"; }

    explicit FieldRefinement(const CFDSim& sim);

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
