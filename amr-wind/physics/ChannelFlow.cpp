#include "amr-wind/physics/ChannelFlow.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/DirectionSelector.H"

namespace amr_wind {
namespace channel_flow {

ChannelFlow::ChannelFlow(CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
{
    
    amrex::ParmParse pp("ChannelFlow");
    pp.query("normal_direction", m_norm_dir);
    
    pp.query("density", m_rho);
    pp.query("re_tau", m_re_tau);
    {
        amrex::Real mu;
        amrex::ParmParse pp("transport");
        pp.query("viscosity", mu);
        // Assumes a boundary layer height of 1.0
        m_utau = mu * m_re_tau / (m_rho * 1.0);
        m_ytau = mu / (m_utau * m_rho);
    }
    pp.query("tke0", m_tke0);
    pp.query("sdr0", m_sdr0);


}

/** Initialize the velocity, density, tke and sdr fields at the beginning of the
 *  simulation.
 */
void ChannelFlow::initialize_fields(
    int level, const amrex::Geometry& geom)
{

    switch (m_norm_dir) {
    case 1:
        initialize_fields(level, geom, YDir(), 1);
        break;
    case 2:
        initialize_fields(level, geom, ZDir(), 2);
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }
}

template <typename IndexSelector>
void ChannelFlow::initialize_fields(
    int level, const amrex::Geometry& geom,
    const IndexSelector& idxOp,
    const int n_idx)
{

    const amrex::Real kappa = m_kappa;
    const amrex::Real y_tau = m_ytau;
    const amrex::Real utau = m_utau;
    auto& velocity = m_repo.get_field("velocity")(level);
    auto& density = m_repo.get_field("density")(level);
    auto& tke = m_repo.get_field("tke")(level);
    auto& sdr = m_repo.get_field("sdr")(level);
    auto& walldist = m_repo.get_field("wall_dist")(level);

    density.setVal(m_rho);
    tke.setVal(m_tke0);
    sdr.setVal(m_sdr0);

    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        const auto& dx = geom.CellSizeArray();
        const auto& problo = geom.ProbLoArray();
        auto vel = velocity.array(mfi);
        auto wd = walldist.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int n_ind = idxOp(i,j,k);
                amrex::Real h = problo[n_idx] + (n_ind + 0.5)*dx[n_idx];
                if (h > 1.0)
                    h = 2.0 - h;
                wd(i,j,k) = h;
                const amrex::Real hp = h/y_tau;
                vel(i,j,k,0) = utau * (1./kappa * std::log(1 + kappa * hp) + 7.8 * ( 1.0 - std::exp(-hp/11.0) - (hp/11.0) * std::exp(-hp/3.0) ) );
                //vel(i,j,k,0) = 22.0;
                vel(i,j,k,1) = 0.0;
                vel(i,j,k,2) = 0.0;
            });
    }
}

void ChannelFlow::post_init_actions() { }

void ChannelFlow::post_advance_work() {

    

    const amrex::Real lam_mu = 1.0e-3;
    auto& den = m_repo.get_field("density");
    auto& tke = m_repo.get_field("tke");
    auto& sdr = m_repo.get_field("sdr");
    const int nlevels = m_repo.num_active_levels();

    // Clip and set values of tke and sdr that are out of bounds
    for (int lev=0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(tke(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& tke_arr = tke(lev).array(mfi);
            const auto& sdr_arr = sdr(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                  const amrex::Real tke_old = tke_arr(i,j,k);
                  const amrex::Real sdr_old = sdr_arr(i,j,k);
                        
                  if ( (tke_old < 0) && (sdr_old < 1.0) ) {
                      tke_arr(i,j,k) = 1e-8;
                      sdr_arr(i,j,k) = rho_arr(i,j,k) * 1e-8 / (1.0*lam_mu);
                  } else if ( (tke_old < 0) ) {
                      tke_arr(i,j,k) = 1.0 * lam_mu * sdr_old / rho_arr(i,j,k);
                  } else if ( (sdr_old < 1.0) ){
                      sdr_arr(i,j,k) = rho_arr(i,j,k) * tke_old / (1.0*lam_mu);
                  }
            });
        }
    }

    amrex::Print() << "Min tke after clip = " << tke(0).min(0) << std::endl;
    amrex::Print() << "Max tke after clip = " << tke(0).max(0) << std::endl;
    amrex::Print() << "Min sdr after clip = " << sdr(0).min(0) << std::endl;
    amrex::Print() << "Max sdr after clip = " << sdr(0).max(0) << std::endl;
    
    tke.fillpatch(m_time.current_time());
    sdr.fillpatch(m_time.current_time());
    
}

} // namespace channel_flow
} // namespace amr_wind
