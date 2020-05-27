#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/diffusion/diffusion.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

ABLWallFunction::ABLWallFunction(const CFDSim& sim)
    : m_sim(sim)
    , m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("normal_direction", m_direction);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_log_law_height);
    } else {
      m_use_fch = true;
    }
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_log_law_height = (geom.ProbLo(m_direction) + 0.5 * geom.CellSize(m_direction));
    } else {

      const auto& geom = m_mesh.Geom();
      amrex::Box const& domain = geom[m_mesh.finestLevel()].Domain();
      const auto dlo = amrex::lbound(domain);
      const auto dhi = amrex::ubound(domain);

      const amrex::Real dz = geom[m_mesh.finestLevel()].CellSize(2);
      m_z_sample_index = dlo.z + std::round((m_log_law_height -
                                             geom[m_mesh.finestLevel()].ProbLo(2))/dz);

      // assuming Z is wall normal direction
      m_ncells_x = dhi.x-dlo.x+1;
      m_ncells_y = dhi.y-dlo.y+1;

      amrex::Real zcellN = geom[m_mesh.finestLevel()].ProbLo(2) + (m_z_sample_index-1)*dz;
      amrex::Real zcellP = geom[m_mesh.finestLevel()].ProbLo(2) + (m_z_sample_index)*dz;

      m_coeff_interp[0] = 1.0 - (m_log_law_height-zcellN)/dz;
      m_coeff_interp[1] = (m_log_law_height-zcellN)/dz;

      amrex::IntVect lo(AMREX_D_DECL(0,0,m_z_sample_index));
      amrex::IntVect hi(AMREX_D_DECL(m_ncells_x-1,m_ncells_y-1,m_z_sample_index));

      m_bx_z_sample.setSmall(lo);
      m_bx_z_sample.setBig(hi);

      m_store_xy_vel.resize(m_bx_z_sample, AMREX_SPACEDIM);

      // amrex::Print() << "Print box" << m_bx_z_sample <<std::endl;
      
    }
}

void ABLWallFunction::update_umean(const FieldPlaneAveraging& pa)
{

  if (m_use_fch) {
    m_umean[m_direction] = 0.0;
    switch (m_direction) {
    case 0:
        m_umean[1] = pa.line_average_interpolated(m_log_law_height, 1);
        m_umean[2] = pa.line_average_interpolated(m_log_law_height, 2);
        break;

    case 1:
        m_umean[0] = pa.line_average_interpolated(m_log_law_height, 0);
        m_umean[2] = pa.line_average_interpolated(m_log_law_height, 2);
        break;

    case 2:
        m_umean[0] = pa.line_average_interpolated(m_log_law_height, 0);
        m_umean[1] = pa.line_average_interpolated(m_log_law_height, 1);
        break;

    default:
        amrex::Abort("Invalid direction specified");
        break;
    }
    m_utau = m_kappa * utils::vec_mag(m_umean.data()) / (
        std::log(m_log_law_height / m_z0));
    amrex::Print() << "Print fch utau" << m_utau <<std::endl;

  } else{
    ComputePlanar();
    m_utau = m_kappa * m_mean_windSpeed/
        (std::log(m_log_law_height/m_z0));
  }
  
}


void ABLWallFunction::ComputePlanar()
{

  amrex::MFItInfo mfi_info{};
  const auto& frepo = m_sim.repo();
  const int nlevels = m_mesh.finestLevel();

  auto geom = m_mesh.Geom();

  auto const& problo = geom[nlevels].ProbLoArray();
  auto const& probhi = geom[nlevels].ProbHiArray();

  const amrex::Real dz = geom[m_mesh.finestLevel()].CellSize(2);
  
  auto& velf = frepo.get_field("velocity", amr_wind::FieldState::New);

  m_store_xy_vel.setVal(0.0);

  auto m_store_xy_vel_arr = m_store_xy_vel.array();

  for (amrex::MFIter mfi(velf(nlevels),mfi_info); mfi.isValid(); ++mfi) {
    const auto& bx = mfi.validbox();
    
    const auto dlo = amrex::lbound(bx);
    const auto dhi = amrex::ubound(bx);

    amrex::Real zminBox = problo[2] + dz*(dlo.z);
    amrex::Real zmaxBox = problo[2] + dz*(dhi.z);

    if((m_log_law_height - zminBox)*(zmaxBox - m_log_law_height)<= 0.0){
      continue;
    }

    auto vel = velf(nlevels).array(mfi);

    const amrex::IntVect lo(dlo.x, dlo.y, m_z_sample_index);
    const amrex::IntVect hi(dhi.x, dhi.y, m_z_sample_index);
    const amrex::Box z_sample_bx(lo, hi);

    amrex::ParallelFor(
          z_sample_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         m_store_xy_vel_arr(i, j, k, 0) = m_coeff_interp[0]*vel(i, j, k-1, 0) + m_coeff_interp[1]*vel(i, j, k, 0);
                         m_store_xy_vel_arr(i, j, k, 1) = m_coeff_interp[0]*vel(i, j, k-1, 1) + m_coeff_interp[1]*vel(i, j, k, 1);
                         m_store_xy_vel_arr(i, j, k, 2) = m_coeff_interp[0]*vel(i, j, k-1, 2) + m_coeff_interp[1]*vel(i, j, k, 2);
              });

  }

  amrex::Real numCells = static_cast<amrex::Real>(m_ncells_x*m_ncells_y);
  
  amrex::ParallelDescriptor::ReduceRealSum(m_store_xy_vel_arr.dataPtr(), numCells*3*sizeof(amrex::Real));

  std::fill(m_umean.begin(), m_umean.end(), 0.0);
  m_mean_windSpeed = 0.0;

  amrex::ParallelFor(
      m_bx_z_sample, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       m_umean[0] += m_store_xy_vel_arr(i, j, k, 0);
                       m_umean[1] += m_store_xy_vel_arr(i, j, k, 1);
                       m_mean_windSpeed += std::sqrt(m_store_xy_vel_arr(i, j, k, 0)*
                                                     m_store_xy_vel_arr(i, j, k, 0) +
                                                     m_store_xy_vel_arr(i, j, k, 1)*
                                                     m_store_xy_vel_arr(i, j, k, 1));
                   });

  m_umean[0] = m_umean[0]/numCells;
  m_umean[1] = m_umean[1]/numCells;
  m_mean_windSpeed = m_mean_windSpeed/numCells;

  amrex::Print() << "Wind speed \t" << m_mean_windSpeed << "\t" << numCells << std::endl;
  
}


ABLVelWallFunc::ABLVelWallFunc(
    Field&, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    diffusion::wall_model_bc(
        velocity, m_wall_func.utau(), utils::vec_mag(m_wall_func.umean().data()),
        rho_state);
}

} // namespace amr_wind
