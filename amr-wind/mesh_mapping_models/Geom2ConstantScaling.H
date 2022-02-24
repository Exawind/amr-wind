#ifndef GEOM2CONSTANTSCALING_H
#define GEOM2CONSTANTSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {
namespace abl_map {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class Geom2ConstantScaling : public MeshMap::Register<Geom2ConstantScaling>
{
public:
    static const std::string identifier() { return "Geom2ConstantScaling"; }

    explicit Geom2ConstantScaling();

    virtual ~Geom2ConstantScaling() = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int, const amrex::Geometry&);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int, const amrex::Geometry&);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&);

private:
    //! User input parameters
    amrex::Vector<amrex::Real> m_sratio{{1.0, 1.0, 1.05}};
    amrex::Vector<amrex::Real> m_delta0{{1.0, 1.0, 1.0}};
    amrex::Vector<amrex::Real> m_transwid{{0.1, 0.1, 0.1}};
    amrex::Vector<amrex::Real> m_transloc{{0.0, 0.0, 10.0}};
    amrex::Vector<int> m_map{{0, 0, 1}};

    amrex::Real m_eps{1e-11};
};

} // namespace abl_map
} // namespace amr_wind

#endif /* Geom2ConstantScaling_H */
