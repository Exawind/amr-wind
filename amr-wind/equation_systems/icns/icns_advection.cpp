#include "amr-wind/equation_systems/icns/icns_advection.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"

#include "AMReX_MacProjector.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace pde {

namespace {

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> bctype,
    amrex::Vector<amrex::Geometry> geom) noexcept
{

    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = amrex::LinOpBCType::Periodic;
        } else {
            auto bc = bctype[amrex::Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_inflow:
            case BC::pressure_outflow: {
                r[dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            default:
                r[dir] = amrex::LinOpBCType::Neumann;
                break;
            };
        }
    }
    return r;
}

} // namespace

MacProjOp::MacProjOp(FieldRepo& repo, bool has_overset, bool variable_density)
    : m_repo(repo)
    , m_options("mac_proj")
    , m_has_overset(has_overset)
    , m_variable_density(variable_density)
{
    amrex::ParmParse pp("incflo");
    pp.query("rho_0", m_rho_0);
}

void MacProjOp::init_projector(const MacProjOp::FaceFabPtrVec& beta) noexcept
{
    m_mac_proj.reset(new amrex::MacProjector(
        m_repo.mesh().Geom(0, m_repo.num_active_levels() - 1)));
    m_mac_proj->initProjector(
        m_options.lpinfo(), beta,
        m_has_overset ? m_repo.get_int_field("mask_cell").vec_const_ptrs()
                      : amrex::Vector<const amrex::iMultiFab*>());

    m_options(*m_mac_proj);

    auto& pressure = m_repo.get_field("p");
    const auto& bctype = pressure.bc_type();

    m_mac_proj->setDomainBC(
        get_projection_bc(
            amrex::Orientation::low, bctype, m_repo.mesh().Geom()),
        get_projection_bc(
            amrex::Orientation::high, bctype, m_repo.mesh().Geom()));

    m_need_init = false;
}

void MacProjOp::init_projector(const amrex::Real beta) noexcept
{
    m_mac_proj.reset(new amrex::MacProjector(
        m_repo.mesh().Geom(0, m_repo.num_active_levels() - 1)));
    m_mac_proj->initProjector(
        m_repo.mesh().boxArray(0, m_repo.num_active_levels() - 1),
        m_repo.mesh().DistributionMap(0, m_repo.num_active_levels() - 1),
        m_options.lpinfo(), beta,
        m_has_overset ? m_repo.get_int_field("mask_cell").vec_const_ptrs()
                      : amrex::Vector<const amrex::iMultiFab*>());

    m_options(*m_mac_proj);

    auto& pressure = m_repo.get_field("p");
    auto& bctype = pressure.bc_type();

    m_mac_proj->setDomainBC(
        get_projection_bc(
            amrex::Orientation::low, bctype, m_repo.mesh().Geom()),
        get_projection_bc(
            amrex::Orientation::high, bctype, m_repo.mesh().Geom()));

    m_need_init = false;
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

void MacProjOp::operator()(const FieldState fstate, const amrex::Real dt)
{
    BL_PROFILE("amr-wind::ICNS::advection_mac_project");
    auto& geom = m_repo.mesh().Geom();
    const auto& pressure = m_repo.get_field("p");
    auto& u_mac = m_repo.get_field("u_mac");
    auto& v_mac = m_repo.get_field("v_mac");
    auto& w_mac = m_repo.get_field("w_mac");
    const auto& density = m_repo.get_field("density", fstate);

    // This will hold (1/rho) on faces
    std::unique_ptr<ScratchField> rho_xf, rho_yf, rho_zf;

    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> rho_face(
        m_repo.num_active_levels());
    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> mac_vec(
        m_repo.num_active_levels());

    // fixme todo clean this up, this was done to replace
    // GetVecOfArrOfConstPtrs() below
    amrex::Vector<amrex::Array<amrex::MultiFab const*, ICNS::ndim>>
        rho_face_const;
    rho_face_const.reserve(m_repo.num_active_levels());

    amrex::Real factor = m_has_overset ? 0.5 * dt : 1.0;

    if (m_variable_density) {

        rho_xf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
        rho_yf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
        rho_zf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);

        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            rho_face[lev][0] = &(*rho_xf)(lev);
            rho_face[lev][1] = &(*rho_yf)(lev);
            rho_face[lev][2] = &(*rho_zf)(lev);

            amrex::average_cellcenter_to_face(
                rho_face[lev], density(lev), geom[lev]);
            for (int idim = 0; idim < ICNS::ndim; ++idim) {
                rho_face[lev][idim]->invert(factor, 0);
            }

            rho_face_const.push_back(GetArrOfConstPtrs(rho_face[lev]));
        }

        if (m_need_init)
            init_projector(rho_face_const);
        else
            m_mac_proj->updateBeta(rho_face_const);

    } else {

        if (m_need_init)
            init_projector(factor / m_rho_0);
        else
            m_mac_proj->updateBeta(factor / m_rho_0);
    }

    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {

        mac_vec[lev][0] = &u_mac(lev);
        mac_vec[lev][1] = &v_mac(lev);
        mac_vec[lev][2] = &w_mac(lev);
    }

    m_mac_proj->setUMAC(mac_vec);

    if (m_has_overset) {
        auto phif = m_repo.create_scratch_field(1, 1, amr_wind::FieldLoc::CELL);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev)
            amrex::average_node_to_cellcenter(
                (*phif)(lev), 0, pressure(lev), 0, 1);

        m_mac_proj->project(
            phif->vec_ptrs(), m_options.rel_tol, m_options.abs_tol);

    } else {
        m_mac_proj->project(m_options.rel_tol, m_options.abs_tol);
    }

    io::print_mlmg_info("MAC_projection", m_mac_proj->getMLMG());
}

} // namespace pde
} // namespace amr_wind
