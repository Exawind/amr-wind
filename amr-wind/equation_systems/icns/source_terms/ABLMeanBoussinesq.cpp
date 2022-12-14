#include "amr-wind/equation_systems/icns/source_terms/ABLMeanBoussinesq.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/wind_energy/ABL.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term for ABL simulations
 *
 *  Reads in the following parameters from `ABLMeanBoussinesq` namespace:
 *
 *  - `reference_temperature` (Mandatory) temperature (`T0`) in Kelvin
 *  - `thermal_expansion_coeff` Optional, default = `1.0 / T0`
 *  - `gravity` acceleration due to gravity (m/s)
 *  - `read_temperature_profile`
 *  - `tprofile_filename`
 */
ABLMeanBoussinesq::ABLMeanBoussinesq(const CFDSim& sim) : m_mesh(sim.mesh())
{

    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_mean_boussinesq_term(this);

    amrex::ParmParse pp_boussinesq_buoyancy("BoussinesqBuoyancy");
    pp_boussinesq_buoyancy.get("reference_temperature", m_ref_theta);

    if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
        pp_boussinesq_buoyancy.get("thermal_expansion_coeff", m_beta);
    } else {
        m_beta = 1.0 / m_ref_theta;
    }

    // FIXME: gravity in `incflo` namespace
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);

    if (pp_boussinesq_buoyancy.contains("read_temperature_profile")) {

        m_const_profile = true;

        std::string tprofile_filename;
        pp_boussinesq_buoyancy.get("tprofile_filename", tprofile_filename);

        read_temperature_profile(tprofile_filename);

    } else {

        mean_temperature_init(abl.abl_statistics().theta_profile());
    }
}

ABLMeanBoussinesq::~ABLMeanBoussinesq() = default;

void ABLMeanBoussinesq::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    const amrex::Real T0 = m_ref_theta;
    const amrex::Real beta = m_beta;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};

    // Mean temperature profile used to compute background forcing term
    //
    // Assumes that the temperature profile is at the cell-centers of the level
    // 0 grid. For finer meshes, it will extrapolate beyond the Level0
    // cell-centers for the lo/hi cells.
    //
    const int idir = m_axis;
    const int nh_max = static_cast<int>(m_theta_ht.size()) - 2;
    const int lp1 = lev + 1;
    const amrex::Real* theights = m_theta_ht.data();
    const amrex::Real* tvals = m_theta_vals.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real temp = T0;
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + (iv[idir] + 0.5) * dx[idir];

        const int il = amrex::min(k / lp1, nh_max);
        const int ir = il + 1;
        temp = tvals[il] +
               ((tvals[ir] - tvals[il]) / (theights[ir] - theights[il])) *
                   (ht - theights[il]);

        const amrex::Real fac = beta * (temp - T0);
        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

void ABLMeanBoussinesq::mean_temperature_init(const FieldPlaneAveraging& tavg)
{
    m_axis = tavg.axis();

    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(tavg.line_centroids().size()));
    m_theta_ht.resize(tavg.line_centroids().size());
    m_theta_vals.resize(tavg.line_average().size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());
    mean_temperature_update(tavg);
}

void ABLMeanBoussinesq::mean_temperature_update(const FieldPlaneAveraging& tavg)
{
    if (m_const_profile) {
        return;
    }
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_average().begin(),
        tavg.line_average().end(), m_theta_vals.begin());
}

void ABLMeanBoussinesq::read_temperature_profile(std::string profile_file_name)
{

    amrex::Vector<amrex::Real> theta_ht, theta_vals;
    std::ifstream infile;
    int n_hts;
    infile.open(profile_file_name.c_str(), std::ios_base::in);
    infile >> n_hts;
    theta_ht.resize(n_hts);
    theta_vals.resize(n_hts);
    m_theta_ht.resize(n_hts);
    m_theta_vals.resize(n_hts);
    for (int i = 0; i < n_hts; i++) {
        infile >> theta_ht[i] >> theta_vals[i];
    }
    infile.close();

    // Now copy to GPU Device memory
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, theta_ht.begin(), theta_ht.end(),
        m_theta_ht.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, theta_vals.begin(), theta_vals.end(),
        m_theta_vals.begin());
}

} // namespace amr_wind::pde::icns
