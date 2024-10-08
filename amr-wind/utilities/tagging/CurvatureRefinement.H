#ifndef CURVATUREREFINEMENT_H
#define CURVATUREREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = CurvatureRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field's curvature value
 *  # tagging/t1.values = 10.0 10.0
 *  ```
 */
class CurvatureRefinement
    : public RefinementCriteria::Register<CurvatureRefinement>
{
public:
    static std::string identifier() { return "CurvatureRefinement"; }

    explicit CurvatureRefinement(const CFDSim& sim);

    ~CurvatureRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_field{nullptr};

    amrex::Vector<amrex::Real> m_curv_value;

    int m_max_lev_field{-1};
};

} // namespace amr_wind

#endif /* CURVATUREREFINEMENT_H */
