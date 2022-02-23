#ifndef ABLSCALING_H
#define ABLSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {
namespace abl_map {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class ABLScaling : public MeshMap::Register<ABLScaling>
{
public:
    static const std::string identifier() { return "ABLScaling"; }

    explicit ABLScaling(const CFDSim& sim);

    virtual ~ABLScaling() = default;

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

    //! User input parameters
    amrex::Vector<amrex::Real> m_sratio{{1.0, 1.0, 1.05}};
    amrex::Vector<amrex::Real> m_delta0{{1.0, 1.0, 1.0}};
    amrex::Vector<int> m_map{{0, 0, 1}};

    amrex::Real m_eps{1e-11};
};

} // namespace abl_map
} // namespace amr_wind

#endif /* ABLSCALING_H */
