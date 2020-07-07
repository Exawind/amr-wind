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
    , m_abl_neutral(true)
{
    amrex::ParmParse pp("ABL");

    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("normal_direction", m_direction);
    pp.queryarr("gravity", m_gravity);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_log_law_height);
    } else {
      m_use_fch = true;
    }

    if (pp.contains("wall_heat_flux"))
    {
      m_abl_neutral = false;
      
      pp.query("wall_heat_flux", m_abl_surface_flux);

      m_heatflux = true;

      if(m_abl_surface_flux<0)
      {
        m_stable = true;
      }else
      {
         m_stable = false;
      }
    } else if(pp.contains("wall_heat_rate"))
    {
      m_abl_neutral = false;
    }
    
}

void ABLWallFunction::init_log_law_height()
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(0);
        m_log_law_height = (geom.ProbLo(m_direction) + 0.5 * geom.CellSize(m_direction));
    } else {

      // m_bx_z_sample.resize(numActiveLevel);
      // m_store_xy_vel.resize(numActiveLevel);
      // m_z_sample_index.resize(numActiveLevel);
      const auto& geom = m_mesh.Geom();

      amrex::Box const& domain = geom[m_mesh.finestLevel()].Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        const amrex::Real dz = geom[m_mesh.finestLevel()].CellSize(2);
        amrex::Real first_cell_height = geom[m_mesh.finestLevel()].ProbLo(2) + 0.5*dz;
        m_z_sample_index = dlo.z +
            static_cast<int> (std::floor((m_log_law_height - first_cell_height)/dz));

        // assuming Z is wall normal direction
        m_ncells_x = dhi.x-dlo.x+1;
        m_ncells_y = dhi.y-dlo.y+1;

        amrex::Real zcellN = first_cell_height + (m_z_sample_index)*dz;
        amrex::Real zcellP = first_cell_height + (m_z_sample_index+1)*dz;

        m_coeff_interp[0] = 1.0 - (m_log_law_height-zcellN)/dz;
        m_coeff_interp[1] = 1.0 - m_coeff_interp[0];

        amrex::Print() << "index of Z sample\t" << m_z_sample_index << "\tInterpolate coeff\t" << m_coeff_interp[0] ;

        amrex::IntVect lo(AMREX_D_DECL(0,0,m_z_sample_index));
        amrex::IntVect hi(AMREX_D_DECL(m_ncells_x-1,m_ncells_y-1,m_z_sample_index));

        m_bx_z_sample.setSmall(lo);
        m_bx_z_sample.setBig(hi);

        // 3 velocity component + potential temperature
        m_store_xy_vel_temp.resize(m_bx_z_sample, 4);

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
  }else{
    computeplanar();
    m_utau = m_kappa * m_mean_windSpeed/
        (std::log(m_log_law_height/m_z0));
    computeusingheatflux();
  }
  
}

void ABLWallFunction::computeplanar()
{

  amrex::MFItInfo mfi_info{};
  const auto& frepo = m_sim.repo();
  const int nlevels = m_mesh.finestLevel();

  auto geom = m_mesh.Geom();

  auto const& problo = geom[nlevels].ProbLoArray();
  auto const& probhi = geom[nlevels].ProbHiArray();

  const amrex::Real dz = geom[m_mesh.finestLevel()].CellSize(2);
  
  auto& velf = frepo.get_field("velocity", amr_wind::FieldState::New);
  auto& tempf = frepo.get_field("temperature", amr_wind::FieldState::New);

  m_store_xy_vel_temp.setVal(0.0);

  auto m_xy_arr = m_store_xy_vel_temp.array();

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
    auto temp = tempf(nlevels).array(mfi);

    const amrex::IntVect lo(dlo.x, dlo.y, m_z_sample_index);
    const amrex::IntVect hi(dhi.x, dhi.y, m_z_sample_index);
    const amrex::Box z_sample_bx(lo, hi);

    amrex::ParallelFor(
          z_sample_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                         m_xy_arr(i, j, k, 0) = m_coeff_interp[0]*vel(i, j, k, 0)
                             + m_coeff_interp[1]*vel(i, j, k+1, 0);
                         m_xy_arr(i, j, k, 1) = m_coeff_interp[0]*vel(i, j, k, 1)
                             + m_coeff_interp[1]*vel(i, j, k+1, 1);
                         m_xy_arr(i, j, k, 2) = m_coeff_interp[0]*vel(i, j, k, 2)
                             + m_coeff_interp[1]*vel(i, j, k+1, 2);
                         m_xy_arr(i, j, k, 3) = m_coeff_interp[0]*temp(i, j, k)
                             + m_coeff_interp[1]*temp(i, j, k+1);

              });

  }

  amrex::Real numCells = static_cast<amrex::Real>(m_ncells_x*m_ncells_y);
  
  amrex::ParallelDescriptor::ReduceRealSum(m_store_xy_vel_temp.dataPtr(),
                                           m_ncells_x*m_ncells_y*3);

  std::fill(m_umean.begin(), m_umean.end(), 0.0);
  m_mean_windSpeed = 0.0;
  m_mean_potTemp = 0.0;

  amrex::ParallelFor(
      m_bx_z_sample, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                          m_umean[0] += m_xy_arr(i, j, k, 0);
                          m_umean[1] += m_xy_arr(i, j, k, 1);
                          m_mean_windSpeed += std::sqrt(m_xy_arr(i, j, k, 0)*
                                                        m_xy_arr(i, j, k, 0) +
                                                        m_xy_arr(i, j, k, 1)*
                                                        m_xy_arr(i, j, k, 1));
                          m_mean_potTemp += m_xy_arr(i, j, k, 3);
                   });

  m_umean[0] = m_umean[0]/numCells;
  m_umean[1] = m_umean[1]/numCells;
  m_mean_windSpeed = m_mean_windSpeed/numCells;
  m_mean_potTemp = m_mean_potTemp/numCells;

}

void ABLWallFunction::computeusingheatflux(){

  if(m_abl_neutral) return;

  amrex::Real g = utils::vec_mag(m_gravity.data());

  amrex::Real zeta = 0.0;
  amrex::Real m_utau_iter = 0.0;
  amrex::Real xzeta = 0.0;
  amrex::Real pi_ = std::acos(-1.0);

  m_psi_m = 0.0;
  m_psi_h = 0.0;

  int iter = 0;
  
  if(m_heatflux){

    if(m_stable) {

      do {
        m_utau_iter = m_utau;
        m_obhukhov_length = -m_utau*m_utau*m_utau*m_mean_potTemp/(m_kappa*g*m_abl_surface_flux);
        zeta = m_log_law_height/m_obhukhov_length;
        m_psi_m = -m_stable_beta*zeta;
        m_psi_h = -m_stable_beta*zeta;
        m_utau = m_kappa * m_mean_windSpeed/(std::log(m_log_law_height/m_z0) - m_psi_m);
        iter +=1;
        
      } while((std::abs(m_utau_iter - m_utau) > 1e-5) && iter <= maxIter);

    } else{

      do {

        m_utau_iter = m_utau;
        m_obhukhov_length = -m_utau*m_utau*m_utau*m_mean_potTemp/(m_kappa*g*m_abl_surface_flux);
        zeta = m_log_law_height/m_obhukhov_length;

        xzeta = std::pow(1.0-m_unstable_gamma*zeta, 0.25);
        
        m_psi_m = 2.0*std::log(0.5*(1.0+xzeta)) + std::log(0.5*(1.0+xzeta*xzeta))
            - 2.0*std::atan(xzeta) + 0.5*pi_;
        m_psi_h = 2.0*std::log(0.5*(1.0+xzeta*xzeta));

        m_utau = m_kappa * m_mean_windSpeed/(std::log(m_log_law_height/m_z0) - m_psi_m);

        iter +=1;

      } while((std::abs(m_utau_iter - m_utau) > 1e-5) && iter <= maxIter);

      amrex::Print()<< " Utau \t" << m_utau << "\t Iter" << iter ;
      
    }

  } else{

    do {
      m_utau_iter = m_utau;
      m_abl_surface_flux = -(m_mean_potTemp-m_abl_surface_temp)*
          m_utau*m_kappa/(std::log(m_log_law_height/m_z0) - m_psi_h);
      m_obhukhov_length = -m_utau*m_utau*m_utau*m_mean_potTemp/(m_kappa*g*m_abl_surface_flux);
    
      zeta = m_log_law_height/m_obhukhov_length;

      if(zeta > 0.0){
        m_psi_m = -m_stable_beta*zeta;
        m_psi_h = -m_stable_beta*zeta;
      } else {
        xzeta = std::pow(1.0-m_unstable_gamma*zeta, 0.25);
        m_psi_m = 2.0*std::log(0.5*(1.0+xzeta)) + std::log(0.5*(1.0+xzeta*xzeta))
            - 2.0*std::atan(xzeta) + 0.5*pi_;
        m_psi_h = 2.0*std::log(0.5*(1.0+xzeta*xzeta));
      }
      m_utau = m_kappa * m_mean_windSpeed/(std::log(m_log_law_height/m_z0) - m_psi_m);
      iter +=1;
    } while((std::abs(m_utau_iter - m_utau) > 1e-5) && iter <= maxIter);
    
  }

  
}


ABLVelWallFunc::ABLVelWallFunc(
    Field&, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    diffusion::wall_model_bc(
        velocity, m_wall_func.utau(), utils::vec_mag(m_wall_func.umean().data()),
        rho_state, m_wall_func.instplanar(), m_wall_func.umean(),
        m_wall_func.meanwindspeed());
}

ABLTempWallFunc::ABLTempWallFunc(
    Field&,
    const ABLWallFunction& wall_fuc)
    : m_wall_func(wall_fuc)
{}

void ABLTempWallFunc::operator()(Field& temperature, const FieldState rho_state){

  diffusion::temp_wall_model_bc(temperature, m_wall_func.heatflux_wall(), rho_state);
}

} // namespace amr_wind
