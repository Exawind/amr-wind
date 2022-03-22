#include <cmath>

#include "amr-wind/mesh_mapping_models/ExpScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace exp_map {

namespace {


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real R(
    const amrex::Real dom_n,
    const amrex::Real delta0,
    const amrex::Real len)
{
    return -1.0*std::log((dom_n)*delta0/len)/(dom_n-1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_fac(
    const amrex::Real ix,
    const amrex::Real x,
    const amrex::Real prob_lo,
    const amrex::Real dom_n,
    const amrex::Real delta0,
    const amrex::Real len,
    const int domap)
{
    // This is the derivative of eval_coord() below
    return (domap == 0)
                ? 1.0
                : (delta0*(ix*R(dom_n,delta0,len)+1.0)*std::exp(R(dom_n,delta0,len)*(ix-1.0))); 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_coord(
    const amrex::Real ix,
    const amrex::Real x,
    const amrex::Real prob_lo,
    const amrex::Real dom_n,
    const amrex::Real delta0,
    const amrex::Real len,
    const int domap)
{
    return (domap == 0)
                ? x
                : (prob_lo + delta0*ix*std::exp(R(dom_n,delta0,len)*(ix-1.0)));
    
    
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_fac_cc(
    const amrex::Real ix,
    const amrex::Real x,
    const amrex::Real prob_lo,
    const amrex::Real dom_n,
    const amrex::Real delta0,
    const amrex::Real len,
    const int domap)
{
    // This is the derivative of eval_coord() below
    return (domap == 0)
                ? 1.0
                : (0.5*delta0*(ix*R(dom_n,delta0,len)+1.0)*std::exp(R(dom_n,delta0,len)*(ix-1.0)))+
                  (0.5*delta0*(ix*R(dom_n,delta0,len)+R(dom_n,delta0,len)+1.0)*std::exp(R(dom_n,delta0,len)*(ix))); 
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_coord_cc(
    const amrex::Real ix,
    const amrex::Real x,
    const amrex::Real prob_lo,
    const amrex::Real dom_n,
    const amrex::Real delta0,
    const amrex::Real len,
    const int domap)
{
    return (domap == 0)
                ? x
                : (prob_lo + 0.5*delta0*ix*std::exp(R(dom_n,delta0,len)*(ix-1.0)))+
                  (prob_lo + 0.5*delta0*(ix+1.0)*std::exp(R(dom_n,delta0,len)*ix));
    
    
}

} // namespace

ExpScaling::ExpScaling()
{
    amrex::ParmParse pp("ExpScaling");
    pp.queryarr("delta0", m_delta0, 0, AMREX_SPACEDIM);
    pp.queryarr("do_map", m_map, 0, AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ExpScaling::create_map(int lev, const amrex::Geometry& geom)
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
void ExpScaling::create_cell_node_map(
    int lev, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> domain_n{
        {len[0]/dx[0], len[1]/dx[1], len[2]/dx[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> x0{
        {prob_lo[0], prob_lo[1], prob_lo[2]}};

    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{
        {m_map[0], m_map[1], m_map[2]}};

    const auto eps = m_eps;

    amrex::Print() << "Map: " << (do_map[0]==1) <<" "<< (do_map[1]==1) <<" "<< (do_map[2]==1)
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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real fac_x = eval_fac_cc(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real fac_y = eval_fac_cc(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real fac_z = eval_fac_cc(iz, z, x0[2] ,domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                scale_fac_cc(i, j, k, 0) = dox ? fac_x : 1.0;
                scale_fac_cc(i, j, k, 1) = doy ? fac_y : 1.0;
                scale_fac_cc(i, j, k, 2) = doz ? fac_z : 1.0;

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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real fac_x = eval_fac(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real fac_y = eval_fac(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real fac_z = eval_fac(iz, z, x0[2] ,domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                scale_fac_nd(i, j, k, 0) = dox ? fac_x : 1.0;
                scale_fac_nd(i, j, k, 1) = doy ? fac_y : 1.0;
                scale_fac_nd(i, j, k, 2) = doz ? fac_z : 1.0;

                scale_detJ_nd(i, j, k) = scale_fac_nd(i, j, k, 0) *
                                         scale_fac_nd(i, j, k, 1) *
                                         scale_fac_nd(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the mesh mapping field on cell faces
 */
void ExpScaling::create_face_map(int lev, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> domain_n{
        {len[0]/dx[0], len[1]/dx[1], len[2]/dx[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> x0{
        {prob_lo[0], prob_lo[1], prob_lo[2]}};
    
    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{
        {m_map[0], m_map[1], m_map[2]}};
    
    const auto eps = m_eps;

    amrex::Print() << "Domain: " << domain_n[0] <<" "<< domain_n[1] <<" "<< domain_n[2]
                       << std::endl;

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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real fac_x = eval_fac(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real fac_y = eval_fac_cc(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real fac_z = eval_fac_cc(iz, z, x0[2] ,domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                scale_fac_xf(i, j, k, 0) = dox ? fac_x : 1.0;
                scale_fac_xf(i, j, k, 1) = doy ? fac_y : 1.0;
                scale_fac_xf(i, j, k, 2) = doz ? fac_z : 1.0;

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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real fac_x = eval_fac_cc(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real fac_y = eval_fac(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real fac_z = eval_fac_cc(iz, z, x0[2] ,domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                scale_fac_yf(i, j, k, 0) = dox ? fac_x : 1.0;
                scale_fac_yf(i, j, k, 1) = doy ? fac_y : 1.0;
                scale_fac_yf(i, j, k, 2) = doz ? fac_z : 1.0;

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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real fac_x = eval_fac_cc(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real fac_y = eval_fac_cc(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real fac_z = eval_fac(iz, z, x0[2] ,domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                scale_fac_zf(i, j, k, 0) = dox ? fac_x : 1.0;
                scale_fac_zf(i, j, k, 1) = doy ? fac_y : 1.0;
                scale_fac_zf(i, j, k, 2) = doz ? fac_z : 1.0;

                scale_detJ_zf(i, j, k) = scale_fac_zf(i, j, k, 0) *
                                         scale_fac_zf(i, j, k, 1) *
                                         scale_fac_zf(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the non-uniform mesh field
 */
void ExpScaling::create_non_uniform_mesh(
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

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> domain_n{
        {len[0]/dx[0], len[1]/dx[1], len[2]/dx[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> delta0{
        {m_delta0[0], m_delta0[1], m_delta0[2]}};

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> x0{
        {prob_lo[0], prob_lo[1], prob_lo[2]}};

    amrex::GpuArray<int, AMREX_SPACEDIM> do_map{
        {m_map[0], m_map[1], m_map[2]}};

    const auto eps = m_eps;

    amrex::Print() << "Domain: " << domain_n[0] <<" "<< domain_n[1] <<" "<< domain_n[2]
                       << std::endl;

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

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real x_non_uni =
                    eval_coord_cc(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real y_non_uni =
                    eval_coord_cc(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real z_non_uni =
                    eval_coord_cc(iz, z, x0[2], domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                nu_coord_cc(i, j, k, 0) = dox ? x_non_uni : x;
                nu_coord_cc(i, j, k, 1) = doy ? y_non_uni : y;
                nu_coord_cc(i, j, k, 2) = doz ? z_non_uni : z;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& nu_coord_nd =
            (*m_non_uniform_coord_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real ix = i;
                amrex::Real iy = j;
                amrex::Real iz = k;

                amrex::Real x_non_uni =
                    eval_coord(ix, x, x0[0], domain_n[0], delta0[0], len[0], do_map[0]);
                amrex::Real y_non_uni =
                    eval_coord(iy, y, x0[1], domain_n[1], delta0[1], len[1], do_map[1]);
                amrex::Real z_non_uni =
                    eval_coord(iz, z, x0[2], domain_n[2], delta0[2], len[2], do_map[2]);

                bool dox = (do_map[0]>0);
                bool doy = (do_map[1]>0);
                bool doz = (do_map[2]>0);

                nu_coord_nd(i, j, k, 0) = dox ? x_non_uni : x;
                nu_coord_nd(i, j, k, 1) = doy ? y_non_uni : y;
                nu_coord_nd(i, j, k, 2) = doz ? z_non_uni : z;
            });
    }
}

  // Used for debugging eval_coord function
  amrex::Real ExpScaling::dump_coord(
                        const amrex::Real ix,
                        const amrex::Real x,
                        const amrex::Real prob_lo,
                        const amrex::Real dom_n,
                        const amrex::Real delta0,
                        const amrex::Real len,
                        const int domap)
  {
    return eval_coord(ix, x, prob_lo, dom_n, delta0, len, domap);
  }

  // Used for debugging eval_fac function
  amrex::Real ExpScaling::dump_fac(
                        const amrex::Real ix,
                        const amrex::Real x,
                        const amrex::Real prob_lo,
                        const amrex::Real dom_n,
                        const amrex::Real delta0,
                        const amrex::Real len,
                        const int domap)
  {
    return eval_fac(ix, x, prob_lo, dom_n, delta0, len, domap);
  }

} // namespace exp_map
} // namespace amr_wind
