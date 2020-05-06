#include <cmath>

#include "ABLFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLFieldInit::ABLFieldInit()
{
    amrex::ParmParse pp("ABL");

    // Temperature variation as a function of height
    pp.getarr("temperature_heights", m_theta_heights);
    pp.getarr("temperature_values", m_theta_values);

    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    
    int num_theta_values = m_theta_heights.size();

    pp.query("perturb_velocity", m_perturb_vel);
    pp.query("perturb_ref_height", m_ref_height);
    pp.query("Uperiods", m_Uperiods);
    pp.query("Vperiods", m_Vperiods);
    pp.query("deltaU", m_deltaU);
    pp.query("deltaV", m_deltaV);

    pp.query("perturb_temperature", m_perturb_theta);

    {
        // TODO: Modify this to accept velocity as a function of height
        // Extract velocity field from incflo
        amrex::ParmParse pp("incflo");
        pp.get("ro_0", m_rho);
        pp.get("ic_u", m_vel[0]);
        pp.get("ic_v", m_vel[1]);
        pp.get("ic_w", m_vel[2]);
    }

    m_thht_d.resize(num_theta_values);
    m_thvv_d.resize(num_theta_values);

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy
#else
    std::memcpy
#endif
    (m_thht_d.data(), m_theta_heights.data(), 
     sizeof(amrex::Real) * num_theta_values);

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy
#else
    std::memcpy
#endif
    (m_thvv_d.data(), m_theta_values.data(), 
     sizeof(amrex::Real) * num_theta_values);
}

void ABLFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& tracer) const
{
    const amrex::Real pi = M_PI;
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    const bool perturb_vel = m_perturb_vel;
    const amrex::Real rho_init = m_rho;
    const amrex::Real umean = m_vel[0];
    const amrex::Real vmean = m_vel[1];
    const amrex::Real wmean = m_vel[2];
    const amrex::Real aval = m_Uperiods * 2.0 * pi / (probhi[1] - problo[1]);
    const amrex::Real bval = m_Vperiods * 2.0 * pi / (probhi[0] - problo[0]);
    const amrex::Real ufac = m_deltaU * std::exp(0.5) / m_ref_height;
    const amrex::Real vfac = m_deltaV * std::exp(0.5) / m_ref_height;
    const amrex::Real ref_height = m_ref_height;

    const int ntvals = m_theta_heights.size();
    const amrex::Real* th = m_thht_d.data();
    const amrex::Real* tv = m_thvv_d.data();

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        density(i, j, k) = rho_init;
        // Mean velocity field
        velocity(i, j, k, 0) = umean;
        velocity(i, j, k, 1) = vmean;
        velocity(i, j, k, 2) = wmean;

        amrex::Real theta = tv[0];
        for (int iz = 0; iz < ntvals - 1; ++iz) {
            if ((z > th[iz]) && (z <= th[iz + 1])) {
                const amrex::Real slope =
                    (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                theta = tv[iz] + (z - th[iz]) * slope;
            }
        }

        // FIXME: Remove first tracer is temperature assumption
        tracer(i, j, k, 0) = theta;

        if (perturb_vel) {
            const amrex::Real xl = x - problo[0];
            const amrex::Real yl = y - problo[1];
            const amrex::Real zl = z / ref_height;
            const amrex::Real damp = std::exp(-0.5 * zl * zl);

            velocity(i, j, k, 0) += ufac * damp * z * std::cos(aval * yl);
            velocity(i, j, k, 1) += vfac * damp * z * std::cos(bval * xl);
        }

        // TODO: Implement temperature perturbations
    });
}

} // namespace amr_wind
