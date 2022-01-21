#ifndef CONSTANTSCALING_H
#define CONSTANTSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {
namespace const_map {

/** Coonstant scaling mesh map
 *  \ingroup mesh_map
 */
class ConstantScaling : public MeshMap::Register<ConstantScaling>
{
public:
    static const std::string identifier() { return "ConstantScaling"; }

    explicit ConstantScaling(const CFDSim& sim);

    virtual ~ConstantScaling() = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell faces
    void create_face_map(int, const amrex::Geometry&) override;

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&) override;

private:
    Field& m_mesh_scale_fac_cc;
    Field& m_mesh_scale_fac_nd;
    Field& m_mesh_scale_fac_xf;
    Field& m_mesh_scale_fac_yf;
    Field& m_mesh_scale_fac_zf;
    Field& m_non_uniform_coord_cc;
    Field& m_non_uniform_coord_nd;

    //! Factor to scale the mesh by
    amrex::Vector<amrex::Real> m_fac{{1.0, 1.0, 1.0}};
};

} // namespace const_map
} // namespace amr_wind

#endif /* CONSTANTSCALING_H */
