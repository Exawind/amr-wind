#include "amr-wind/core/MeshMap.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {

void MeshMap::declare_mapping_fields(const CFDSim& sim, int nghost)
{

    // declare nodal, cell-centered, and face-centered mesh mapping array
    m_mesh_scale_fac_cc = &(sim.repo().declare_cc_field(
        "mesh_scaling_factor_cc", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_nd = &(sim.repo().declare_nd_field(
        "mesh_scaling_factor_nd", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_xf = &(sim.repo().declare_xf_field(
        "mesh_scaling_factor_xf", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_yf = &(sim.repo().declare_yf_field(
        "mesh_scaling_factor_yf", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_zf = &(sim.repo().declare_zf_field(
        "mesh_scaling_factor_zf", AMREX_SPACEDIM, nghost, 1));

    // declare nodal, cell-centered, and face-centered mesh mapping detJ array
    m_mesh_scale_detJ_cc =
        &(sim.repo().declare_cc_field("mesh_scaling_detJ_cc", 1, nghost, 1));
    m_mesh_scale_detJ_nd =
        &(sim.repo().declare_nd_field("mesh_scaling_detJ_nd", 1, nghost, 1));
    m_mesh_scale_detJ_xf =
        &(sim.repo().declare_xf_field("mesh_scaling_detJ_xf", 1, nghost, 1));
    m_mesh_scale_detJ_yf =
        &(sim.repo().declare_yf_field("mesh_scaling_detJ_yf", 1, nghost, 1));
    m_mesh_scale_detJ_zf =
        &(sim.repo().declare_zf_field("mesh_scaling_detJ_zf", 1, nghost, 1));

    // declare nodal and cell-centered non-uniform mesh
    m_non_uniform_coord_cc = &(sim.repo().declare_cc_field(
        "non_uniform_coord_cc", AMREX_SPACEDIM, nghost, 1));
    m_non_uniform_coord_nd = &(sim.repo().declare_nd_field(
        "non_uniform_coord_nd", AMREX_SPACEDIM, nghost, 1));

    // TODO: Create BCNoOP fill patch operators for mesh scaling fields ?
}

void MeshMap::setup_interp_arrays(int lev, const amrex::Geometry& geom) 
{
    const auto& prob_lo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    amrex::Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    const int Nx = dhi.x - dlo.x + 1;
    const int Ny = dhi.y - dlo.y + 1;
    const int Nz = dhi.z - dlo.z + 1;

    amrex::Vector<amrex::Real> x_uni;
    amrex::Vector<amrex::Real> y_uni;
    amrex::Vector<amrex::Real> z_uni;
    x_uni.resize(Nx, 0.0);
    y_uni.resize(Ny, 0.0);
    z_uni.resize(Nz, 0.0);

    amrex::Vector<amrex::Real> x_nu;
    amrex::Vector<amrex::Real> y_nu;
    amrex::Vector<amrex::Real> z_nu;
    x_nu.resize(Nx, 0.0);
    y_nu.resize(Ny, 0.0);
    z_nu.resize(Nz, 0.0);

    amr_wind::Field const* nu_coord_cc = m_non_uniform_coord_cc;

    for (amrex::MFIter mfi((*m_mesh_scale_fac_cc)(lev)); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real const> nu_cc =
            ((*nu_coord_cc)(lev).array(mfi));
	// Create the uniform vector coordinates
	for (int i = 0; i < Nx; i++) {
	  x_uni[i] = prob_lo[0] + (i + 0.5) * dx[0];
	}
	for (int j = 0; j < Ny; j++) {
	  y_uni[j] = prob_lo[1] + (j + 0.5) * dx[1];
	}
	for (int k = 0; k<Nz; k++) {
	  z_uni[k] = prob_lo[2] + (k + 0.5) * dx[2];
	}
	amrex::Loop(bx, [=, &x_nu, &y_nu, &z_nu](int i, int j, int k) noexcept {
	    if ((j == 0) && (k == 0) && (i >= 0) && (i < Nx)){
	      x_nu[i] = nu_cc(i,j,k,0);
	    }
	    if ((i == 0) && (k == 0) && (j >= 0) && (j < Ny)){
	      y_nu[j] = nu_cc(i,j,k,1);
	    }
	    if ((i == 0) && (j == 0) && (k >= 0) && (k < Nz)){
	      z_nu[k] = nu_cc(i,j,k,2);
	    }
	  });
    }

    // Copy over to storage vectors
    m_x_uni = x_uni;
    m_y_uni = y_uni;
    m_z_uni = z_uni;

    m_x_nu = x_nu;
    m_y_nu = y_nu;
    m_z_nu = z_nu;
}

amrex::Real MeshMap::interp_nonunif_to_unif(amrex::Real x_nonunif, int dir)
{
    namespace interp = ::amr_wind::interp;
    amrex::Real ans(0.0);
    switch(dir) {
    case 0: 
        ans = interp::linear(m_x_nu, m_x_uni, x_nonunif);
        break;
    case 1: 
        ans = interp::linear(m_y_nu, m_y_uni, x_nonunif);
        break;
    case 2: 
        ans = interp::linear(m_z_nu, m_z_uni, x_nonunif);
        break;
    }
    return ans;
}
  
} // namespace amr_wind
