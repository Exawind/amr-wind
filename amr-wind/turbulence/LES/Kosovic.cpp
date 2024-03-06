#include <cmath>

#include "amr-wind/turbulence/LES/Kosovic.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/nonLinearSum.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/divergence.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
Kosovic<Transport>::Kosovic(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
    , m_Nij(sim.repo().declare_field("Nij", 9, 1, 1))
    , m_divNij(sim.repo().declare_field("divNij", 3, 1, 1))
{
    amrex::ParmParse pp("Kosovic");
    pp.query("Cb", m_Cb);
    m_Cs = std::sqrt(8 * (1 + m_Cb) / (27 * 3.14 * 3.14));
    m_C1 = std::sqrt(960) * m_Cb / (7 * (1 + m_Cb) * m_Sk);
    m_C2 = m_C1;
    amrex::Print() << "C1:" << m_C1 << std::endl;
    pp.query("refMOL", m_refMOL);
    if (std::abs(m_refMOL) > 1000)
        amrex::Print()
            << " Length scale assumes Neutral Stratification with Ref MOL:"
            << m_refMOL << std::endl;
    else if (m_refMOL < 0)
        amrex::Print()
            << " Length scale suggests Unstable  Stratification with Ref MOL:"
            << m_refMOL << std::endl;
    else
        amrex::Print()
            << " Length scale suggests Stable  Stratification with Ref MOL::"
            << m_refMOL << std::endl;
    pp.query("surfaceRANS", m_surfaceRANS);
    if (m_surfaceRANS) {
        m_surfaceFactor = 1;
        pp.query("switchLoc", m_switchLoc);
        pp.query("surfaceRANSExp", m_surfaceRANSExp);
        amrex::Print() << " Enabling RANS model in the surface layer upto "
                       << m_switchLoc << "  m" << std::endl;
        amrex::Print() << " Exponential term power coefficient:"
                       << m_surfaceRANSExp << std::endl;
    } else {
        m_surfaceFactor = 0;
        amrex::Print() << " No RANS model in the surface layer" << std::endl;
    }
    pp.query("writeTerms", m_writeTerms);
    if (m_writeTerms) {
        this->m_sim.io_manager().register_io_var("Nij");
        this->m_sim.io_manager().register_io_var("divNij");
    } else
        amrex::Print() << " Not writing the non-linear terms" << std::endl;
    pp.query("LESOff", m_LESTurnOff);
    amrex::Print() << " Turning off SGS model above " << m_LESTurnOff << "   m"
                   << std::endl;
}

template <typename Transport>
void Kosovic<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;

    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);
    // Non-linear component Nij is computed here and goes into Body Forcing
    // fvm::gradient(m_gradU,vel);
    fvm::nonlinearsum(m_Nij, vel);
    fvm::divergence(m_divNij, m_Nij);
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);
            const auto& divNijLevel = (this->m_divNij)(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    const amrex::Real x3 = problo[2] + (k + 0.5) * dz;
                    // if(i==76 && j==76)
                    // amrex::Print()<<i<<"   "<<j<<"  "<<k<<"  "<<x3<<"
                    // "<<dz<<std::endl;
                    const amrex::Real fmu = std::exp(-x3 / m_switchLoc);
                    amrex::Real phiM = 1;
                    if (m_refMOL < 0)
                        phiM = std::pow(1 - 16 * x3 / m_refMOL, -0.25);
                    else
                        phiM = 1 + 5 * x3 / m_refMOL;
                    const amrex::Real ransL =
                        std::pow(0.41 * (k + 1) * dz / phiM, 2);
                    amrex::Real turnOff = std::exp(-x3 / m_LESTurnOff);
                    // if(x3>=m_LESTurnOff)
                    //   turnOff=0;
                    // const amrex::Real
                    // masonScaleInv=1.0/(std::pow(0.41*(k+1)*dz,2))+1/smag_factor;
                    // amrex::Real
                    // hybridScale=amrex::min(smag_factor,std::pow(0.41*(k+1)*dz,2));
                    // amrex::Real
                    // lscale=m_surfaceFactor*hybridScale+(1-m_surfaceFactor)*smag_factor;
                    amrex::Real lscale =
                        m_surfaceFactor *
                            (std::pow(1 - fmu, m_surfaceRANSExp) * smag_factor +
                             std::pow(fmu, m_surfaceRANSExp) * ransL) +
                        (1 - m_surfaceFactor) * smag_factor;
                    mu_arr(i, j, k) *= rho * lscale * turnOff;
                    // lscale=m_surfaceFactor*hybridScale+(1-m_surfaceFactor)*smag_factor*0.25*m_C1;
                    lscale = m_surfaceFactor *
                                 (std::pow(1 - fmu, m_surfaceRANSExp) *
                                      smag_factor * 0.25 * m_C1 +
                                  std::pow(fmu, m_surfaceRANSExp) * ransL) +
                             (1 - m_surfaceFactor) * smag_factor * 0.25 * m_C1;
                    divNijLevel(i, j, k, 0) *= rho * lscale * turnOff;
                    divNijLevel(i, j, k, 1) *= rho * lscale * turnOff;
                    divNijLevel(i, j, k, 2) *= rho * lscale * turnOff;
                });
        }
    }

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void Kosovic<Transport>::update_alphaeff(Field& alphaeff)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    //    const auto& geom_vec = repo.mesh().Geom();
    amrex::Real muCoeff = 3;
    if (m_refMOL > 0) muCoeff = 1;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& muturb_arr = mu_turb(lev).array(mfi);
            const auto& alphaeff_arr = alphaeff(lev).array(mfi);
            const auto& lam_diff_arr = (*lam_alpha)(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    alphaeff_arr(i, j, k) =
                        lam_diff_arr(i, j, k) + muCoeff * muturb_arr(i, j, k);
                });
        }
    }
    alphaeff.fillpatch(this->m_sim.time().current_time());
}
template <typename Transport>
void Kosovic<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType Kosovic<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cs", this->m_Cs}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Kosovic);

} // namespace amr_wind
