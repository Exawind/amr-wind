#ifndef EXPSCALING_H
#define EXPSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {
namespace exp_map {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class ExpScaling : public MeshMap::Register<ExpScaling>
{
public:
    static const std::string identifier() { return "ExpScaling"; }

    explicit ExpScaling();

    virtual ~ExpScaling() = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int, const amrex::Geometry&);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int, const amrex::Geometry&);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&);

        //! Used for debugging eval_coord function
    amrex::Real dump_coord(
			   const int ,
			   const amrex::Real ,
			   const amrex::Real ,
			   const amrex::Real ,
			   const amrex::Real ,
			   const amrex::Real ,
			   const int);
    //! Used for debugging eval_fac function
    amrex::Real dump_fac(
			 const int ,
			 const amrex::Real ,
			 const amrex::Real ,
			 const amrex::Real ,
			 const amrex::Real ,
			 const amrex::Real ,
			 const int);

private:
    //! User input parameters
    amrex::Vector<amrex::Real> m_delta0{{1.0, 1.0, 1.0}};
    amrex::Vector<int> m_map{{0, 0, 1}};

    amrex::Real m_eps{1e-11};
};

} // namespace abl_map
} // namespace amr_wind

#endif /* ExpScaling_H */
