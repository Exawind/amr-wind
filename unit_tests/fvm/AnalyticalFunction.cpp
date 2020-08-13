#include "AnalyticalFunction.H"

namespace amr_wind_tests {

AnalyticalFunction::AnalyticalFunction(int degree, int ncomponents)
    : m_degree(degree), m_ncomponents(ncomponents)
{
    // fixme do we want random variables assigned here?
    // fixme do we want different coeffs for each component?
    m_coeffs.assign(std::pow(degree + 1, AMREX_SPACEDIM), 2.0);
}

AnalyticalFunction::~AnalyticalFunction() = default;

amrex::Real AnalyticalFunction::phi_eval(
    amrex::Real x, amrex::Real y, amrex::Real z, int component) const
{
    amrex::Vector<amrex::Real> phi;
    phi.assign(m_ncomponents, 0.0);

    int ind = 0.0;
    for (int i = 0; i <= m_degree; ++i) {
        for (int j = 0; j <= m_degree; ++j) {
            for (int k = 0; k <= m_degree; ++k) {
                for (int n = 0; n < m_ncomponents; ++n) {
                    phi[n] += m_coeffs[ind] * std::pow(x, i) * std::pow(y, j) *
                              std::pow(z, k);
                }
                ++ind;
            }
        }
    }

    return phi[component];
}

amrex::Real AnalyticalFunction::dphi_eval(
    amrex::Real x, amrex::Real y, amrex::Real z, int component) const
{
    amrex::Vector<amrex::Real> dphi;
    dphi.assign(AMREX_SPACEDIM * m_ncomponents, 0.0);

    int ind = 0;
    for (int i = 0; i <= m_degree; ++i) {
        for (int j = 0; j <= m_degree; ++j) {
            for (int k = 0; k <= m_degree; ++k) {
                for (int n = 0; n < m_ncomponents; ++n) {
                    dphi[AMREX_SPACEDIM * n + 0] +=
                        i * m_coeffs[ind] * std::pow(x, amrex::max(i - 1, 0)) *
                        std::pow(y, j) * std::pow(z, k);
                    dphi[AMREX_SPACEDIM * n + 1] +=
                        j * m_coeffs[ind] * std::pow(x, i) *
                        std::pow(y, amrex::max(j - 1, 0)) * std::pow(z, k);
                    dphi[AMREX_SPACEDIM * n + 2] +=
                        k * m_coeffs[ind] * std::pow(x, i) * std::pow(y, j) *
                        std::pow(z, amrex::max(k - 1, 0));
                }
                ++ind;
            }
        }
    }

    return dphi[component];
}

amrex::Real
AnalyticalFunction::laplacian(amrex::Real x, amrex::Real y, amrex::Real z) const
{
    amrex::Real lap = 0.0;

    int ind = 0;
    for (int i = 0; i <= m_degree; ++i) {
        for (int j = 0; j <= m_degree; ++j) {
            for (int k = 0; k <= m_degree; ++k) {
                lap += i * amrex::max(i - 1, 0) * m_coeffs[ind] *
                       std::pow(x, amrex::max(i - 2, 0)) * std::pow(y, j) *
                       std::pow(z, k);
                lap += j * amrex::max(j - 1, 0) * m_coeffs[ind] *
                       std::pow(x, i) * std::pow(y, amrex::max(j - 2, 0)) *
                       std::pow(z, k);
                lap += k * amrex::max(k - 1, 0) * m_coeffs[ind] *
                       std::pow(x, i) * std::pow(y, j) *
                       std::pow(z, amrex::max(k - 2, 0));

                ++ind;
            }
        }
    }

    return lap;
}

amrex::Real AnalyticalFunction::divergence(
    amrex::Real x, amrex::Real y, amrex::Real z) const
{
    AMREX_ALWAYS_ASSERT(m_ncomponents == AMREX_SPACEDIM);

    return dphi_eval(x, y, z, 0) + dphi_eval(x, y, z, 4) +
           dphi_eval(x, y, z, 8);
}

amrex::Real AnalyticalFunction::strainrate(
    amrex::Real x, amrex::Real y, amrex::Real z) const
{
    AMREX_ALWAYS_ASSERT(m_ncomponents == AMREX_SPACEDIM);

    const int ncomp = AMREX_SPACEDIM * AMREX_SPACEDIM;
    amrex::Vector<amrex::Real> gradvec(ncomp);

    for (int n = 0; n < ncomp; ++n) gradvec[n] = dphi_eval(x, y, z, n);

    return std::sqrt(
        2.0 * gradvec[0] * gradvec[0] + 2.0 * gradvec[4] * gradvec[4] +
        2.0 * gradvec[8] * gradvec[8] +
        (gradvec[1] + gradvec[3]) * (gradvec[1] + gradvec[3]) +
        (gradvec[5] + gradvec[7]) * (gradvec[5] + gradvec[7]) +
        (gradvec[6] + gradvec[2]) * (gradvec[6] + gradvec[2]));
}

} // namespace amr_wind_tests
