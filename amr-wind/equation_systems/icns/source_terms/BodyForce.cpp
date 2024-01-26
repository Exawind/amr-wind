#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <cstddef>
#include <ios>

namespace amr_wind::pde::icns {

/** Body Force
 */
BodyForce::BodyForce(const CFDSim& sim) : m_time(sim.time()), m_mesh(sim.mesh())
{

    // Read the geostrophic wind speed vector (in m/s)
    amrex::ParmParse pp("BodyForce");
    pp.query("type", m_type);
    m_type = amrex::toLower(m_type);

    if (m_type == "height-varying") {
        pp.get("bodyforce-file", m_bforce_file);
        read_bforce_profile(m_bforce_file);
    } else {
        pp.getarr("magnitude", m_body_force);
        if (m_type == "oscillatory") {
            pp.get("angular_frequency", m_omega);
        }
    }
}

BodyForce::~BodyForce() = default;

void BodyForce::read_bforce_profile(const std::string& filename)
{

    std::ifstream infile;
    amrex::Vector<amrex::Real> bforce_x;
    amrex::Vector<amrex::Real> bforce_y;
    amrex::Vector<amrex::Real> bforce_hts;

    infile.open(filename.c_str(), std::ios_base::in);
    infile >> m_bforce_profile_nhts;

    bforce_x.resize(m_bforce_profile_nhts);
    bforce_y.resize(m_bforce_profile_nhts);
    bforce_hts.resize(m_bforce_profile_nhts);

    m_prof_x.resize(m_bforce_profile_nhts);
    m_prof_y.resize(m_bforce_profile_nhts);
    m_ht.resize(m_bforce_profile_nhts);
    for (size_t i = 0; i < m_bforce_profile_nhts; i++) {
        infile >> bforce_hts[i] >> bforce_x[i] >> bforce_y[i];
    }
    infile.close();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_x.begin(), bforce_x.end(),
        m_prof_x.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_y.begin(), bforce_y.end(),
        m_prof_y.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_hts.begin(), bforce_hts.end(),
        m_ht.begin());
}

void BodyForce::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& time = m_time.current_time();
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    const int lp1 = lev + 1;
    const int nh_max = (int)m_prof_x.size() - 2;

    const amrex::Real* force_ht = m_ht.data();
    const amrex::Real* force_x = m_prof_x.data();
    const amrex::Real* force_y = m_prof_y.data();

    if (m_type == "height-varying") {

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::IntVect iv(i, j, k);
                const amrex::Real ht = problo[2] + (iv[2] + 0.5) * dx[2];
                const int il = amrex::min(k / lp1, nh_max);
                const int ir = il + 1;
                amrex::Real fx;
                amrex::Real fy;

                fx = force_x[il] + ((force_x[ir] - force_x[il]) /
                                    (force_ht[ir] - force_ht[il])) *
                                       (ht - force_ht[il]);

                fy = force_y[il] + ((force_y[ir] - force_y[il]) /
                                    (force_ht[ir] - force_ht[il])) *
                                       (ht - force_ht[il]);

                src_term(i, j, k, 0) += fx;
                src_term(i, j, k, 1) += fy;
            });

    } else {

        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
            {m_body_force[0], m_body_force[1], m_body_force[2]}};

        amrex::Real coeff =
            (m_type == "oscillatory") ? std::cos(m_omega * time) : 1.0;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                src_term(i, j, k, 0) += coeff * forcing[0];
                src_term(i, j, k, 1) += coeff * forcing[1];
                src_term(i, j, k, 2) += coeff * forcing[2];
            });
    }
}

} // namespace amr_wind::pde::icns
