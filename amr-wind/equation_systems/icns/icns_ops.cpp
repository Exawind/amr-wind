#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "AMReX_MacProjector.H"
#include "amr-wind/convection/mac_projection.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"

namespace amr_wind {
namespace pde {

//
// Computes the following decomposition:
//
//    u + c*grad(phi)/ro = u*  with  div(ep*u) = 0
//
// Inputs:
//
//   u_mac,v_mac,w_mac = the MAC velocity field to be projected
//   density           = the cell-centered density
//
// Outputs:
//
//  u_mac,v_mac,w_mac = the PROJECTED MAC velocity field
//
// Notes:
//
//  phi, the projection auxiliary function, is computed by solving
//
//       div(ep*grad(phi)/rho) = div(ep * u*)
//

void advection_mac_project(FieldRepo& repo, const FieldState fstate, const bool has_overset)
{
    BL_PROFILE("amr-wind::ICNS::advection_mac_project");
    auto& geom = repo.mesh().Geom();
    auto& u_mac = repo.get_field("u_mac");
    auto& v_mac = repo.get_field("v_mac");
    auto& w_mac = repo.get_field("w_mac");
    auto& density = repo.get_field("density", fstate);

    // This will hold (1/rho) on faces
    auto rho_xf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
    auto rho_yf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
    auto rho_zf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);

    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> rho_face(
        repo.num_active_levels());
    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> mac_vec(
        repo.num_active_levels());

    // fixme todo clean this up, this was done to replace
    // GetVecOfArrOfConstPtrs() below
    amrex::Vector<amrex::Array<amrex::MultiFab const*, ICNS::ndim>>
        rho_face_const;
    rho_face_const.reserve(repo.num_active_levels());

    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        rho_face[lev][0] = &(*rho_xf)(lev);
        rho_face[lev][1] = &(*rho_yf)(lev);
        rho_face[lev][2] = &(*rho_zf)(lev);

        amrex::average_cellcenter_to_face(
            rho_face[lev], density(lev), geom[lev]);
        for (int idim = 0; idim < ICNS::ndim; ++idim) {
            rho_face[lev][idim]->invert(1.0, 0);
        }

        rho_face_const.push_back(GetArrOfConstPtrs(rho_face[lev]));

        mac_vec[lev][0] = &u_mac(lev);
        mac_vec[lev][1] = &v_mac(lev);
        mac_vec[lev][2] = &w_mac(lev);
    }

    amr_wind::MLMGOptions options("mac_proj");

    // If we want to set max_coarsening_level we have to send it in to the
    // constructor
    amrex::LPInfo lp_info;
    lp_info.setMaxCoarseningLevel(options.max_coarsen_level);

    // Perform MAC projection
    amrex::MacProjector macproj(
        mac_vec, amrex::MLMG::Location::FaceCenter,
        rho_face_const, amrex::MLMG::Location::FaceCenter,
        amrex::MLMG::Location::CellCenter,
        repo.mesh().Geom(0, repo.num_active_levels() - 1), lp_info, {},
        amrex::MLMG::Location::CellCenter,
        (has_overset ? repo.get_int_field("mask_cell").vec_const_ptrs()
                     : amrex::Vector<const amrex::iMultiFab*>()));

    // get the bc types from pressure field
    auto& bctype = repo.get_field("p").bc_type();

    macproj.setDomainBC(
        mac::get_projection_bc(
            amrex::Orientation::low, bctype, repo.mesh().Geom()),
        mac::get_projection_bc(
            amrex::Orientation::high, bctype, repo.mesh().Geom()));

    macproj.project(options.rel_tol, options.abs_tol);

    io::print_mlmg_info("MAC_projection", macproj.getMLMG());
}

} // namespace pde
} // namespace amr_wind
