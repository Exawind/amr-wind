#ifndef VORTICITYREFINEMENT_H
#define VORTICITYREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
class Field;

/** AMR refinement using vorticity
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = v1
 *  tagging/v1.type = VorticityRefinement
 *  # Tagging based on value
 *  tagging/v1.values = 10.0 10.0

 *  ```
 */
class VorticityRefinement
    : public RefinementCriteria::Register<VorticityRefinement>
{
public:
    static std::string identifier() { return "VorticityRefinement"; }

    VorticityRefinement(const CFDSim& sim);

    virtual ~VorticityRefinement() = default;

    //! Read input file and initialize boxarray used to refine each level
    virtual void initialize(const std::string& key) override;

    virtual void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_vel;

    amrex::Vector<amrex::Real> m_vort_value;

    int m_max_lev_field;
};

} // namespace amr_wind

#endif /* VORTICITYREFINEMENT_H */
