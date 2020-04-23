#include <cmath>

#include "Smagorinsky.H"
#include "TurbModelDefs.H"
#include "derive_K.H"

#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template<typename Transport>
Smagorinsky<Transport>::Smagorinsky(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", m_Cs);
}

template <typename Transport>
void Smagorinsky<Transport>::update_turbulent_viscosity(const FieldState fstate)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_turbulent_viscosity")

    auto& mu_turb = this->mu_turb();
    auto& repo = mu_turb.repo();
    auto& vel = repo.get_field("velocity", fstate);
    auto& den = repo.get_field("density", fstate);
    auto& geom_vec = repo.mesh().Geom();
    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;

    const int nlevels = repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& domain = geom.Domain();

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;

        const amrex::Real idx = 1.0 / dx;
        const amrex::Real idy = 1.0 / dy;
        const amrex::Real idz = 1.0 / dz;

        for (amrex::MFIter mfi(mu_turb(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox(mu_turb.num_grow());
            const auto& mu_arr = mu_turb(lev).array(mfi);
            const auto& vel_arr = vel(lev).const_array(mfi);
            const auto& rho_arr = den(lev).const_array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real rho = rho_arr(i, j, k);
                    const amrex::Real sr = incflo_strainrate<StencilInterior>(
                        i, j, k, idx, idy, idz, vel_arr);
                    mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                });

            // TODO: Check if the following is correct for `foextrap` BC types
            const auto& bxi = mfi.tilebox();
            int idim = 0;
            if (!geom.isPeriodic(idim)) {
                if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                    amrex::IntVect low(bxi.smallEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxlo = amrex::Box(low, hi).grow({0, 1, 1});

                    amrex::ParallelFor(
                        bxlo, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilILO>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }

                if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                    amrex::IntVect low(bxi.bigEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxhi = amrex::Box(low, hi).grow({0, 1, 1});

                    amrex::ParallelFor(
                        bxhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilIHI>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }
            } // if (!geom.isPeriodic)

            idim = 1;
            if (!geom.isPeriodic(idim)) {
                if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                    amrex::IntVect low(bxi.smallEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxlo = amrex::Box(low, hi).grow({1, 0, 1});

                    amrex::ParallelFor(
                        bxlo, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilJLO>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }

                if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                    amrex::IntVect low(bxi.bigEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxhi = amrex::Box(low, hi).grow({1, 0, 1});

                    amrex::ParallelFor(
                        bxhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilJHI>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }
            } // if (!geom.isPeriodic)

            idim = 2;
            if (!geom.isPeriodic(idim)) {
                if (bxi.smallEnd(idim) == domain.smallEnd(idim)) {
                    amrex::IntVect low(bxi.smallEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxlo = amrex::Box(low, hi).grow({1, 1, 0});

                    amrex::ParallelFor(
                        bxlo, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilKLO>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }

                if (bxi.bigEnd(idim) == domain.bigEnd(idim)) {
                    amrex::IntVect low(bxi.bigEnd());
                    amrex::IntVect hi(bxi.bigEnd());
                    int sm = low[idim];
                    low.setVal(idim, sm);
                    hi.setVal(idim, sm);

                    auto bxhi = amrex::Box(low, hi).grow({1, 1, 0});

                    amrex::ParallelFor(
                        bxhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            const amrex::Real rho = rho_arr(i, j, k);
                            const amrex::Real sr = incflo_strainrate<StencilKHI>(
                                i, j, k, idx, idy, idz, vel_arr);
                            mu_arr(i, j, k) = rho * Cs_sqr * ds_sqr * sr;
                        });
                }
            } // if (!geom.isPeriodic)
        }
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Smagorinsky);

} // namespace amr_wind
