#ifndef MESHMAP_H
#define MESHMAP_H

#include "amr-wind/core/CollMgr.H"
#include "amr-wind/core/Factory.H"
#include "amr-wind/core/Field.H"

#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

namespace amr_wind {

class CFDSim;

/**
 *  \defgroup mesh_map Mesh mapping models
 *
 *  AMR-Wind representation of different mesh mapping models
 *
 *  In AMR-Wind, different mesh mappings are implemented using the MeshMap
 * class.
 */

/** Abstract representation of different physics
 *  \ingroup physics
 *
 *  This class defines an abstract API that represents the notion of some
 *  mesh mapping that will be used to scale the mesh. The most common use-case
 *  for this class is to perform RANS simulations.
 */
class MeshMap : public Factory<MeshMap>
{
public:
    static const std::string base_identifier() { return "MeshMap"; }

    virtual ~MeshMap() = default;

    //! declare mesh mapping fields
    void declare_mapping_fields(const CFDSim&, int);

    //! Construct mesh scaling field
    virtual void create_map(int, const amrex::Geometry&) = 0;

    //! Construct mesh scaling field on cell centers and nodes
    virtual void create_cell_node_map(int, const amrex::Geometry&) = 0;

    //! Construct mesh scaling field on cell faces
    virtual void create_face_map(int, const amrex::Geometry&) = 0;

    //! Construct the non-uniform mesh field
    virtual void create_non_uniform_mesh(int, const amrex::Geometry&) = 0;

protected:
    Field* m_mesh_scale_fac_cc{nullptr};
    Field* m_mesh_scale_fac_nd{nullptr};
    Field* m_mesh_scale_fac_xf{nullptr};
    Field* m_mesh_scale_fac_yf{nullptr};
    Field* m_mesh_scale_fac_zf{nullptr};
    Field* m_non_uniform_coord_cc{nullptr};
    Field* m_non_uniform_coord_nd{nullptr};
};

} // namespace amr_wind

#endif /* MESHMAP_H */
