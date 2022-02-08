#ifndef CONSTANTMAP_H
#define CONSTANTMAP_H

#include "amr-wind/core/MeshMap.H"

namespace amr_wind {
namespace const_map {

/** Coonstant scaling mesh map
 *  \ingroup mesh_map
 */
class ConstantMap : public MeshMap::Register<ConstantMap>
{
public:
    static const std::string identifier() { return "ConstantMap"; }

    explicit ConstantMap();

    virtual ~ConstantMap() = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&);

private:
    //! Factor to scale the mesh by
    amrex::Vector<amrex::Real> m_fac{{1.0, 1.0, 1.0}};
};

} // namespace const_map
} // namespace amr_wind

#endif /* CONSTANTMAP_H */
