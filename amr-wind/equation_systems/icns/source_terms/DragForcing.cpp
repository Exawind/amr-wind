#include "amr-wind/equation_systems/icns/source_terms/DragForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind::pde::icns {

DragForcing::DragForcing(const CFDSim& sim)
  :m_sim(sim),
   m_mesh(sim.mesh()),
   m_velocity(sim.repo().get_field("velocity")),
   m_terrainBlank(sim.repo().get_field("terrainBlank")),
   m_terrainDrag(sim.repo().get_field("terrainDrag"))
   
{
    const auto& abl = m_sim.physics_manager().get<amr_wind::ABL>();
    const VelPlaneAveraging& fa_velocity=abl.abl_statistics().vel_profile_coarse();
    gpu_vel_ht.resize(fa_velocity.line_centroids().size());
    gpu_vel_vals.resize(fa_velocity.line_average().size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, fa_velocity.line_centroids().begin(),
                     fa_velocity.line_centroids().end(), gpu_vel_ht.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, fa_velocity.line_average().begin(),
                     fa_velocity.line_average().end(), gpu_vel_vals.begin());
    /*for(int ii=0;ii<gpu_vel_ht.size();++ii)
      amrex::Print()<<"H:"<<gpu_vel_ht[ii]
                    <<" u:"<<gpu_vel_vals[3*ii]
                    <<" v:"<<gpu_vel_vals[3*ii+1]
                    <<" w:"<<gpu_vel_vals[3*ii+2]
                    <<std::endl;*/
}
  
DragForcing::~DragForcing() = default;

void DragForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
      m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const auto blank = m_terrainBlank(lev).const_array(mfi);
    const auto drag = m_terrainDrag(lev).const_array(mfi);
    const auto& geom_vec = m_mesh.Geom();
    const auto& geom = geom_vec[lev];   
    const auto& dx = geom.CellSize();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    //const amrex::Real dy = geom.CellSize()[1];
    //const amrex::Real dz = geom.CellSize()[2];
    //Real xi = (x - xhi_sponge_start) / (ProbHiArr[0] - xhi_sponge_start);
    //rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
    const amrex::Real gpu_drag=m_drag;
    const amrex::Real gpu_spongeStrength=m_spongeStrength;
    const amrex::Real gpu_spongeDensity=m_spongeDensity;
    const amrex::Real gpu_startX=(1-m_spongePercentX/100.0)*prob_hi[0];
    const amrex::Real gpu_startY=(1-m_spongePercentY/100.0)*prob_hi[1];
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real ux=vel(i,j,k,0);
      const amrex::Real uy=vel(i,j,k,1);
      const amrex::Real uz=vel(i,j,k,2);
      const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
      const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
      amrex::Real xdamping=0;
      amrex::Real ydamping=0;
      if(x1>gpu_startX){
	amrex::Real xi=(x1 - gpu_startX) / (prob_hi[0] - gpu_startX);
	xdamping=gpu_spongeStrength * xi * xi ;
      }
        if(x2>gpu_startY){
        amrex::Real yi=(x2 - gpu_startY) / (prob_hi[1] - gpu_startY);
	ydamping=gpu_spongeStrength * yi * yi;
      }
      const amrex::Real m=std::sqrt(ux*ux+uy*uy+uz*uz);
      amrex::Real Cd=gpu_drag/dx[0];
      amrex::Real gpu_spongeVelX=0.0;
      amrex::Real gpu_spongeVelY=0.0;
      amrex::Real gpu_spongeVelZ=0.0;
      amrex::Real residual=1000;
      amrex::Real height_error=0.0;
      for(int ii=0;ii<gpu_vel_ht.size();++ii)
	{
	  height_error=std::abs(x3-gpu_vel_ht[ii]);
	  if(height_error<residual)
	    {
	      residual=height_error;
	      gpu_spongeVelX=gpu_vel_vals[3*ii];
	      gpu_spongeVelY=gpu_vel_vals[3*ii+1];
	      gpu_spongeVelZ=gpu_vel_vals[3*ii+2];
	    }
	}
      // Terrain Drag
      amrex::Real kappa=0.41;
      amrex::Real ustar=std::sqrt(ux*ux+uy*uy)*kappa/std::log((x3+0.1)/0.1);
      amrex::Real Dxz= -ustar*ustar*ux/(1e-5+std::sqrt(ux*ux+uy*uy))/dx[2];
      amrex::Real Dyz= -ustar*ustar*uy/(1e-5+std::sqrt(ux*ux+uy*uy))/dx[2];
      // Adjusting Cd for momentum
      amrex::Real CdM=std::min(Cd*5.0/(m+1e-5),100.0);
      src_term(i, j, k, 0) -= (CdM*m*ux*blank(i,j,k)+Dxz*drag(i,j,k)+(xdamping+ydamping)*(ux-gpu_spongeDensity*gpu_spongeVelX));
      src_term(i, j, k, 1) -= (CdM*m*uy*blank(i,j,k)+Dyz*drag(i,j,k)+(xdamping+ydamping)*(uy-gpu_spongeDensity*gpu_spongeVelY));
      src_term(i, j, k, 2) -= (CdM*m*uz*blank(i,j,k)+(xdamping+ydamping)*(uz-gpu_spongeDensity*gpu_spongeVelZ));
    });
}

} // namespace amr_wind::pde::icns
