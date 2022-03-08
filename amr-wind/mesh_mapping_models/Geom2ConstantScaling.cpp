#include <cmath>

#include "amr-wind/mesh_mapping_models/Geom2ConstantScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace geom2constant_map {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real x_geom(
    const amrex::Real x,
    const amrex::Real sratio,
    const amrex::Real delta0)
{
    return delta0/sratio/(1.0-sratio)*(1.0-std::pow(sratio, x+1.0))-delta0/sratio;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real x_const(
    const amrex::Real x,
    const amrex::Real slope,
    const amrex::Real xinv,
    const amrex::Real lmatch)
{
    return slope*(x - xinv) + lmatch;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real x_inv(
    const amrex::Real trans_loc,
    const amrex::Real sratio,
    const amrex::Real delta0)
{
    return int(std::log(1.0-(trans_loc-delta0/sratio)*sratio*(1.0-sratio)/delta0)/std::log(sratio));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_fac(
    const amrex::Real x,
    const amrex::Real trans_loc,
    const amrex::Real trans_width,
    const amrex::Real sratio,
    const amrex::Real delta0,
    const amrex::Real dxscale,
    const int domap)
{
    // This is the derivative of eval_coord() below
    const amrex::Real xinv = x_inv(trans_loc, sratio, delta0);
    const amrex::Real h = 0.5*(1.0+std::tanh((x*dxscale-xinv)/trans_width));
    const amrex::Real hprime = 0.5*std::pow(std::cosh((x*dxscale-xinv)/trans_width),-2.0)/trans_width;
    const amrex::Real lmatch = x_geom(xinv, sratio, delta0);
    const amrex::Real slope = lmatch-x_geom(xinv-1.0, sratio, delta0);
    const amrex::Real xprime_geom = delta0/sratio/(sratio-1.0)*std::pow(sratio, x-1.0)*std::log(sratio);

    if (domap==0) return 1.0;
    return ((1.0-h)*xprime_geom - hprime*x_geom(x,sratio,delta0) + h*slope + hprime*x_const(x,slope,xinv,lmatch))*dxscale;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_coord(
    const amrex::Real x,
    const amrex::Real trans_loc,
    const amrex::Real trans_width,
    const amrex::Real sratio,
    const amrex::Real delta0,
    const amrex::Real dxscale,
    const int domap)
{
    const amrex::Real xinv = x_inv(trans_loc, sratio, delta0);
    const amrex::Real h = 0.5*(1.0+std::tanh((x*dxscale-xinv)/trans_width));
    const amrex::Real lmatch = x_geom(xinv, sratio, delta0);
    const amrex::Real slope = lmatch-x_geom(xinv-1.0, sratio, delta0);

    if (domap==0) return x;
    return x_geom(x*dxscale,sratio,delta0)*(1.0-h) + h*x_const(x*dxscale,slope,xinv,lmatch);
}

} // namespace

Geom2ConstantScaling::Geom2ConstantScaling()
{
    amrex::ParmParse pp("Geom2ConstantScaling");
    pp.queryarr("sratio", m_sratio, 0, AMREX_SPACEDIM);
    pp.queryarr("delta0", m_delta0, 0, AMREX_SPACEDIM);
    pp.queryarr("transwidth", m_transwid, 0, AMREX_SPACEDIM);
    pp.queryarr("translocation", m_transloc, 0, AMREX_SPACEDIM);
    pp.queryarr("do_map", m_map, 0, AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void Geom2ConstantScaling::create_map(int lev, const amrex::Geometry& geom)
{
    create_cell_node_map(lev, geom);
    create_face_map(lev, geom);
    create_non_uniform_mesh(lev, geom);
    if (lev==0) {
        setup_interp_arrays(lev, geom);
    }

}

/** Construct the mesh mapping field on cell centers and nodes
 */
void Geom2ConstantScaling::create_cell_node_map(
    int lev, const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sratio{
        {m_sratio[0], m_sratio[1], m_sratio[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transwid{
        {m_transwid[0], m_transwid[1], m_transwid[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transloc{
        {m_transloc[0], m_transloc[1], m_transloc[2]}};
    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{{m_map[0], m_map[1], m_map[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxscale{{1.0, 1.0, 1.0}};
    for (int i=0; i<AMREX_SPACEDIM; i++) {
      if (do_map[i]) {
	dxscale[i] = 1.0/dx[i];
      }
    }

    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    amrex::Print() << "Transition Width:    " << transwid[0] <<" "<< transwid[1] <<" "<< transwid[2]
                       << std::endl;

    amrex::Print() << "Transition Location:    " << transloc[0] <<" "<< transloc[1] <<" "<< transloc[2]
                       << std::endl;

    for (amrex::MFIter mfi((*m_mesh_scale_fac_cc)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc =
            (*m_mesh_scale_fac_cc)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_cc =
            (*m_mesh_scale_detJ_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
	      amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
	      amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = 
    		    eval_fac(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real fac_y = 
                    eval_fac(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real fac_z = 
                    eval_fac(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real fac_x_below = do_map[0] ? delta0[0]*dxscale[0] : 1.0;
		const amrex::Real fac_y_below = do_map[1] ? delta0[1]*dxscale[1] : 1.0;
		const amrex::Real fac_z_below = do_map[2] ? delta0[2]*dxscale[2] : 1.0;
                scale_fac_cc(i, j, k, 0) = x_below ? fac_x_below : fac_x;
                scale_fac_cc(i, j, k, 1) = y_below ? fac_y_below : fac_y;
                scale_fac_cc(i, j, k, 2) = z_below ? fac_z_below : fac_z;

                scale_detJ_cc(i, j, k) = scale_fac_cc(i, j, k, 0) *
                                         scale_fac_cc(i, j, k, 1) *
                                         scale_fac_cc(i, j, k, 2);
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd =
            (*m_mesh_scale_fac_nd)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_nd =
            (*m_mesh_scale_detJ_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      amrex::Real x = prob_lo[0] + i * dx[0];
	      amrex::Real y = prob_lo[1] + j * dx[1];
	      amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x = 
                    eval_fac(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real fac_y = 
                    eval_fac(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real fac_z = 
                    eval_fac(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real fac_x_below = do_map[0] ? delta0[0]*dxscale[0] : 1.0;
		const amrex::Real fac_y_below = do_map[1] ? delta0[1]*dxscale[1] : 1.0;
		const amrex::Real fac_z_below = do_map[2] ? delta0[2]*dxscale[2] : 1.0;
                scale_fac_nd(i, j, k, 0) = x_below ? fac_x_below : fac_x;
                scale_fac_nd(i, j, k, 1) = y_below ? fac_y_below : fac_y;
                scale_fac_nd(i, j, k, 2) = z_below ? fac_z_below : fac_z;

                scale_detJ_nd(i, j, k) = scale_fac_nd(i, j, k, 0) *
                                         scale_fac_nd(i, j, k, 1) *
                                         scale_fac_nd(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the mesh mapping field on cell faces
 */
void Geom2ConstantScaling::create_face_map(int lev, const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sratio{
        {m_sratio[0], m_sratio[1], m_sratio[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transwid{
        {m_transwid[0], m_transwid[1], m_transwid[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transloc{
        {m_transloc[0], m_transloc[1], m_transloc[2]}};
    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{{m_map[0], m_map[1], m_map[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxscale{{1.0, 1.0, 1.0}};
    for (int i=0; i<AMREX_SPACEDIM; i++) {
      if (do_map[i]) {
	dxscale[i] = 1.0/dx[i];
      }
    }

    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    for (amrex::MFIter mfi((*m_mesh_scale_fac_xf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_xf =
            (*m_mesh_scale_fac_xf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_xf =
            (*m_mesh_scale_detJ_xf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      amrex::Real x = prob_lo[0] + i * dx[0];
	      amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
	      amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = 
                    eval_fac(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real fac_y = 
                    eval_fac(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real fac_z = 
                    eval_fac(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real fac_x_below = do_map[0] ? delta0[0]*dxscale[0] : 1.0;
		const amrex::Real fac_y_below = do_map[1] ? delta0[1]*dxscale[1] : 1.0;
		const amrex::Real fac_z_below = do_map[2] ? delta0[2]*dxscale[2] : 1.0;
                scale_fac_xf(i, j, k, 0) = x_below ? fac_x_below : fac_x;
                scale_fac_xf(i, j, k, 1) = y_below ? fac_y_below : fac_y;
                scale_fac_xf(i, j, k, 2) = z_below ? fac_z_below : fac_z;

                scale_detJ_xf(i, j, k) = scale_fac_xf(i, j, k, 0) *
                                         scale_fac_xf(i, j, k, 1) *
                                         scale_fac_xf(i, j, k, 2);
            });
    }

    for (amrex::MFIter mfi((*m_mesh_scale_fac_yf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_yf =
            (*m_mesh_scale_fac_yf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_yf =
            (*m_mesh_scale_detJ_yf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	      amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
	      amrex::Real y = prob_lo[1] + j * dx[1];
	      amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = 
                    eval_fac(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real fac_y = 
                    eval_fac(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real fac_z = 
                    eval_fac(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real fac_x_below = do_map[0] ? delta0[0]*dxscale[0] : 1.0;
		const amrex::Real fac_y_below = do_map[1] ? delta0[1]*dxscale[1] : 1.0;
		const amrex::Real fac_z_below = do_map[2] ? delta0[2]*dxscale[2] : 1.0;
                scale_fac_yf(i, j, k, 0) = x_below ? fac_x_below : fac_x;
                scale_fac_yf(i, j, k, 1) = y_below ? fac_y_below : fac_y;
                scale_fac_yf(i, j, k, 2) = z_below ? fac_z_below : fac_z;

                scale_detJ_yf(i, j, k) = scale_fac_yf(i, j, k, 0) *
                                         scale_fac_yf(i, j, k, 1) *
                                         scale_fac_yf(i, j, k, 2);
            });
    }

    for (amrex::MFIter mfi((*m_mesh_scale_fac_zf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_zf =
            (*m_mesh_scale_fac_zf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_zf =
            (*m_mesh_scale_detJ_zf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	        amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
	        amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
	        amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x = 
                    eval_fac(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real fac_y = 
                    eval_fac(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real fac_z = 
                    eval_fac(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real fac_x_below = do_map[0] ? delta0[0]*dxscale[0] : 1.0;
		const amrex::Real fac_y_below = do_map[1] ? delta0[1]*dxscale[1] : 1.0;
		const amrex::Real fac_z_below = do_map[2] ? delta0[2]*dxscale[2] : 1.0;
                scale_fac_zf(i, j, k, 0) = x_below ? fac_x_below : fac_x;
                scale_fac_zf(i, j, k, 1) = y_below ? fac_y_below : fac_y;
                scale_fac_zf(i, j, k, 2) = z_below ? fac_z_below : fac_z;

                scale_detJ_zf(i, j, k) = scale_fac_zf(i, j, k, 0) *
                                         scale_fac_zf(i, j, k, 1) *
                                         scale_fac_zf(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the non-uniform mesh field
 */
void Geom2ConstantScaling::create_non_uniform_mesh(
    int lev, const amrex::Geometry& geom)
{
    amrex::Vector<amrex::Real> probhi_physical{{0.0, 0.0, 0.0}};
    {
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_physical")) {
            pp.getarr("prob_hi_physical", probhi_physical);
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                probhi_physical[d] = geom.ProbHiArray()[d];
            }
        }
    }

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sratio{
        {m_sratio[0], m_sratio[1], m_sratio[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transwid{
        {m_transwid[0], m_transwid[1], m_transwid[2]}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> transloc{
        {m_transloc[0], m_transloc[1], m_transloc[2]}};
    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{{m_map[0], m_map[1], m_map[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxscale{{1.0, 1.0, 1.0}};
    for (int i=0; i<AMREX_SPACEDIM; i++) {
      if (do_map[i]) {
	dxscale[i] = 1.0/dx[i];
      }
    }

    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {probhi_physical[0] - prob_lo[0], probhi_physical[1] - prob_lo[1],
         probhi_physical[2] - prob_lo[2]}};

    for (amrex::MFIter mfi((*m_non_uniform_coord_cc)(lev)); mfi.isValid();
         ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& nu_coord_cc =
            (*m_non_uniform_coord_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
	        amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
	        amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
	        amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real x_non_uni =
                    eval_coord(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real y_non_uni =
                    eval_coord(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real z_non_uni =
                    eval_coord(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real cc_x_below = do_map[0] ? prob_lo[0] + (i + 0.5) * delta0[0] : x;
		const amrex::Real cc_y_below = do_map[1] ? prob_lo[1] + (j + 0.5) * delta0[1] : y;
		const amrex::Real cc_z_below = do_map[2] ? prob_lo[2] + (k + 0.5) * delta0[2] : z;
                nu_coord_cc(i, j, k, 0) = x_below ? cc_x_below : x_non_uni;
                nu_coord_cc(i, j, k, 1) = y_below ? cc_y_below : y_non_uni;
                nu_coord_cc(i, j, k, 2) = z_below ? cc_z_below : z_non_uni;

            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& nu_coord_nd =
            (*m_non_uniform_coord_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
  	        amrex::Real x = prob_lo[0] + i * dx[0];
		amrex::Real y = prob_lo[1] + j * dx[1];
		amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real x_non_uni =
                    eval_coord(x, transloc[0], transwid[0], sratio[0], delta0[0], dxscale[0], do_map[0]);
                amrex::Real y_non_uni =
                    eval_coord(y, transloc[1], transwid[1], sratio[1], delta0[1], dxscale[1], do_map[1]);
                amrex::Real z_non_uni =
                    eval_coord(z, transloc[2], transwid[2], sratio[2], delta0[2], dxscale[2], do_map[2]);

		bool x_below = x < prob_lo[0];
		bool y_below = y < prob_lo[1];
		bool z_below = z < prob_lo[2];
		const amrex::Real nd_x_below = do_map[0] ? prob_lo[0] + i * delta0[0] : x;
		const amrex::Real nd_y_below = do_map[1] ? prob_lo[1] + j * delta0[1] : y;
		const amrex::Real nd_z_below = do_map[2] ? prob_lo[2] + k * delta0[2] : z;

                nu_coord_nd(i, j, k, 0) = x_below ? nd_x_below : x_non_uni;
                nu_coord_nd(i, j, k, 1) = y_below ? nd_y_below : y_non_uni;
                nu_coord_nd(i, j, k, 2) = z_below ? nd_z_below : z_non_uni;
            });
    }
}

  // Used for debugging eval_coord function
  amrex::Real Geom2ConstantScaling::dump_coord(
					       const amrex::Real x,
					       const amrex::Real trans_loc,
					       const amrex::Real trans_width,
					       const amrex::Real sratio,
					       const amrex::Real delta0,
					       const amrex::Real len,
					       const int domap)
  {
    return eval_coord(x, trans_loc, trans_width, sratio, delta0, len, domap);
  }

  // Used for debugging eval_fac function
  amrex::Real Geom2ConstantScaling::dump_fac(
					     const amrex::Real x,
					     const amrex::Real trans_loc,
					     const amrex::Real trans_width,
					     const amrex::Real sratio,
					     const amrex::Real delta0,
					     const amrex::Real len,
					     const int domap)
  {
    return eval_fac(x, trans_loc, trans_width, sratio, delta0, len, domap);
  }

} // namespace geom2constant_map
} // namespace amr_wind
