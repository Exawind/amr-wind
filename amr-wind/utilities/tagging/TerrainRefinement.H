#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
// class CFDSim;
class Field;
// class IntField;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = TerrainRefinement
 *  tagging/t1.level = 3
 *  tagging/t1.vertical_distance = 1´500
 *  ```
 */
class TerrainRefinement : public RefinementCriteria::Register<TerrainRefinement>
{
public:
    static std::string identifier() { return "TerrainRefinement"; }

    explicit TerrainRefinement(const CFDSim& sim);

    ~TerrainRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void operator()(
        const int level,
        amrex::TagBoxArray& tags,
        const amrex::Real time,
        const int ngrow) override;

private:
    const CFDSim& m_sim;

    Field* m_terrain_height{nullptr};
    IntField* m_terrain_blank{nullptr};
    int m_max_lev{0};
    amrex::Real m_vertical_distance{0.0};
    amrex::RealBox m_tagging_box;
};

} // namespace amr_wind

#endif /* TERRAINREFINEMENT_H */
