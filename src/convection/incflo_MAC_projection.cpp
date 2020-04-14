#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>
#include <AMReX_Vector.H>

#include <mac_projection.H>
#include <MLMGOptions.H>

using namespace amrex;

Array<amrex::LinOpBCType,AMREX_SPACEDIM>
mac::get_projection_bc (Orientation::Side side, GpuArray<BC, AMREX_SPACEDIM*2> bctype, Vector<Geometry> geom) noexcept
{

    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = bctype[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            case BC::mass_inflow:
            case BC::slip_wall:
            case BC::no_slip_wall:
            case BC::wall_model:
            {
                r[dir] = LinOpBCType::Neumann;
                break;
            }
            default:
                amrex::Abort("mac get_projection_bc: undefined BC type");
            };
        }
    }
    return r;
}

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
void 
mac::apply_MAC_projection (amr_wind::FieldRepo& repo,
                           amr_wind::FieldState fstate)
{
    BL_PROFILE("incflo::apply_MAC_projection()");

    auto& u_mac = repo.get_field("u_mac");
    auto& v_mac = repo.get_field("v_mac");
    auto& w_mac = repo.get_field("w_mac");
    auto& density = repo.get_field("density", fstate);
    auto& geom = repo.mesh().Geom();

//    if (m_verbose > 2) amrex::Print() << "MAC Projection:\n";

    // This will hold (1/rho) on faces
    auto rho_xf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
    auto rho_yf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
    auto rho_zf = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);

    Vector<Array<MultiFab*,AMREX_SPACEDIM>> rho_face(repo.num_active_levels());
    Vector<Array<MultiFab*,AMREX_SPACEDIM>> mac_vec(repo.num_active_levels());

    //fixme todo clean this up, this was done to replace GetVecOfArrOfConstPtrs() below
    Vector<Array<MultiFab const*,AMREX_SPACEDIM>> rho_face_const;
    rho_face_const.reserve(repo.num_active_levels());

    for (int lev=0; lev < repo.num_active_levels(); ++lev)
    {
        rho_face[lev][0] = &(*rho_xf)(lev);
        rho_face[lev][1] = &(*rho_yf)(lev);
        rho_face[lev][2] = &(*rho_zf)(lev);

        amrex::average_cellcenter_to_face(rho_face[lev], density(lev), geom[lev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rho_face[lev][idim]->invert(1.0, 0);
        }

        rho_face_const.push_back(GetArrOfConstPtrs(rho_face[lev]));

        mac_vec[lev][0] = &u_mac(lev);
        mac_vec[lev][1] = &v_mac(lev);
        mac_vec[lev][2] = &w_mac(lev);
    }

    amr_wind::MLMGOptions options("mac_proj");
    //
    // If we want to set max_coarsening_level we have to send it in to the constructor
    //
    LPInfo lp_info;
    lp_info.setMaxCoarseningLevel(options.max_coarsen_level);

    //
    // Perform MAC projection
    //
    MacProjector macproj(mac_vec,
                         rho_face_const,
                         repo.mesh().Geom(0,repo.num_active_levels()-1),
                         lp_info);

    // get the bc types from pressure field
    auto& bctype = repo.get_field("p").bc_type();
    
    macproj.setDomainBC(get_projection_bc(Orientation::low, bctype, repo.mesh().Geom()),
                        get_projection_bc(Orientation::high, bctype, repo.mesh().Geom()));

    macproj.project(
        options.rel_tol, options.abs_tol, MLMG::Location::FaceCentroid);
}
