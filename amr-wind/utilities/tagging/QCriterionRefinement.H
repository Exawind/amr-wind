#ifndef QCRITERIONREFINEMENT_H
#define QCRITERIONREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

/** AMR refinement using Q criterion
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = qc1
 *  tagging/qc1.type = QCriterionRefinement
 *  # Tagging based on value
 *  tagging/qc1.values = 10.0 10.0

 *  ```
 */
class QCriterionRefinement
    : public RefinementCriteria::Register<QCriterionRefinement>
{
public:
    static std::string identifier() { return "QCriterionRefinement"; }

    explicit QCriterionRefinement(const CFDSim& sim);

    virtual ~QCriterionRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_vel;

    amrex::Vector<amrex::Real> m_qc_value;

    int m_max_lev_field;

    bool m_nondim{true};
};

} // namespace amr_wind

#endif /* QCRITERIONREFINEMENT_H */
