#ifndef GRADIENTMAGREFINEMENT_H
#define GRADIENTMAGREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = GradientMagRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field's gradient magnitude value
 *  # tagging/t1.values = 10.0 10.0
<<<<<<< HEAD
=======

>>>>>>> updated comments for Curvature and Gradient magnitude based refinement criteria
 */
class GradientMagRefinement
    : public RefinementCriteria::Register<GradientMagRefinement>
{
public:
    static std::string identifier() { return "GradientMagRefinement"; }

    GradientMagRefinement(const CFDSim& sim);

    virtual ~GradientMagRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_field;

    amrex::Vector<amrex::Real> m_gradmag_value;

    int m_max_lev_field;
};

} // namespace amr_wind

#endif /* GRADIENTMAGREFINEMENT_H */
