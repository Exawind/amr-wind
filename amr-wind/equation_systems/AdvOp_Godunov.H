#ifndef ADVOP_GODUNOV_H
#define ADVOP_GODUNOV_H

#include <type_traits>

#include "amr-wind/convection/Godunov.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "hydro_utils.H"

namespace amr_wind::pde {

/** Godunov advection operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::Godunov,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
{
    AdvectionOp(
        PDEFields& fields_in,
        bool /* has_overset */,
        bool /* variable density */,
        bool /* mesh mapping */)
        : fields(fields_in)
        , density(fields_in.repo.get_field("density"))
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {
        amrex::ParmParse pp("incflo");
        pp.query("godunov_type", godunov_type);
        pp.query("godunov_use_forces_in_trans", godunov_use_forces_in_trans);
        if (pp.contains("use_ppm") || pp.contains("use_limiter")) {
            amrex::Abort(
                "Godunov: use_ppm and use_limiter are deprecated. Please "
                "update input file");
        }

        if (amrex::toLower(godunov_type) == "plm") {
            godunov_scheme = godunov::scheme::PLM;
        } else if (amrex::toLower(godunov_type) == "ppm") {
            godunov_scheme = godunov::scheme::PPM;
        } else if (amrex::toLower(godunov_type) == "ppm_nolim") {
            godunov_scheme = godunov::scheme::PPM_NOLIM;
        } else if (amrex::toLower(godunov_type) == "bds") {
            godunov_scheme = godunov::scheme::BDS;
            advection_type = "BDS";
        } else if (
            amrex::toLower(godunov_type) == "weno" ||
            amrex::toLower(godunov_type) == "weno_js") {
            godunov_scheme = godunov::scheme::WENOJS;
        } else if (amrex::toLower(godunov_type) == "weno_z") {
            godunov_scheme = godunov::scheme::WENOZ;
        } else {
            amrex::Print()
                << "For godunov_type select between plm, ppm, "
                   "ppm_nolim, bds, weno_js, and weno_z: it defaults to ppm"
                << std::endl;
            godunov_scheme = godunov::scheme::PPM;
        }
        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(PDE::ndim, 1);
    }

    void preadvect(const FieldState /*unused*/, const amrex::Real /*unused*/) {}

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        static_assert(
            PDE::ndim == 1, "Invalid number of components for scalar");
        auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        // cppcheck-suppress constVariable
        auto& conv_term = fields.conv_term;
        const auto& dof_field = fields.field.state(fstate);

        auto flux_x =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::ZFACE);
        auto face_x =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto face_y =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto face_z =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::ZFACE);

        // only needed if multiplying by rho below
        const auto& den = density.state(fstate);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                auto rho_arr = den(lev).array(mfi);
                auto tra_arr = dof_field(lev).array(mfi);
                amrex::FArrayBox rhotracfab;
                amrex::Array4<amrex::Real> rhotrac;

                if (PDE::multiply_rho) {
                    auto rhotrac_box =
                        amrex::grow(bx, fvm::Godunov::nghost_state);
                    rhotracfab.resize(rhotrac_box, PDE::ndim);
                    rhotrac = rhotracfab.array();

                    amrex::ParallelFor(
                        rhotrac_box, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(
                            int i, int j, int k, int n) noexcept {
                            rhotrac(i, j, k, n) =
                                rho_arr(i, j, k) * tra_arr(i, j, k, n);
                        });
                }

                if ((godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                    (godunov_scheme == godunov::scheme::WENOJS) ||
                    (godunov_scheme == godunov::scheme::WENOZ)) {
                    amrex::FArrayBox tmpfab(amrex::grow(bx, 1), PDE::ndim * 14);

                    godunov::compute_fluxes(
                        lev, bx, PDE::ndim, (*flux_x)(lev).array(mfi),
                        (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                        (PDE::multiply_rho ? rhotrac : tra_arr),
                        u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi),
                        src_term(lev).const_array(mfi),
                        dof_field.bcrec_device().data(), iconserv.data(),
                        tmpfab.dataPtr(), geom, dt, godunov_scheme);
                } else if (
                    (godunov_scheme == godunov::scheme::PPM) ||
                    (godunov_scheme == godunov::scheme::PLM) ||
                    (godunov_scheme == godunov::scheme::BDS)) {
                    amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 1);
                    tmpfab.setVal<amrex::RunOn::Device>(0.0);
                    const auto& divu = tmpfab.array();
                    const bool is_velocity = false;
                    const bool known_edge_state = false;
                    const bool godunov_use_ppm =
                        godunov_scheme == godunov::scheme::PPM;
                    const bool fluxes_are_area_weighted = false;
                    HydroUtils::ComputeFluxesOnBoxFromState(
                        bx, PDE::ndim, mfi,
                        (PDE::multiply_rho ? rhotrac : tra_arr),
                        AMREX_D_DECL(
                            (*flux_x)(lev).array(mfi),
                            (*flux_y)(lev).array(mfi),
                            (*flux_z)(lev).array(mfi)),
                        AMREX_D_DECL(
                            (*face_x)(lev).array(mfi),
                            (*face_y)(lev).array(mfi),
                            (*face_z)(lev).array(mfi)),
                        known_edge_state,
                        AMREX_D_DECL(
                            u_mac(lev).const_array(mfi),
                            v_mac(lev).const_array(mfi),
                            w_mac(lev).const_array(mfi)),
                        divu, src_term(lev).const_array(mfi), geom[lev], dt,
                        dof_field.bcrec(), dof_field.bcrec_device().data(),
                        iconserv.data(), godunov_use_ppm,
                        godunov_use_forces_in_trans, is_velocity,
                        fluxes_are_area_weighted, advection_type);
                } else {
                    amrex::Abort("Invalid godunov scheme");
                }
                amrex::Gpu::streamSynchronize();
            }
        }

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
            repo.num_active_levels());
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            fluxes[lev][0] = &(*flux_x)(lev);
            fluxes[lev][1] = &(*flux_y)(lev);
            fluxes[lev][2] = &(*flux_z)(lev);
        }

        // In order to enforce conservation across coarse-fine boundaries we
        // must be sure to average down the fluxes before we use them
        for (int lev = repo.num_active_levels() - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1], rr,
                geom[lev - 1]);
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();

                godunov::compute_advection(
                    lev, bx, PDE::ndim, conv_term(lev).array(mfi),
                    (*flux_x)(lev).array(mfi), (*flux_y)(lev).array(mfi),
                    (*flux_z)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    iconserv.data(), geom);
            }
        }
    }

    PDEFields& fields;
    Field& density;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    amrex::Gpu::DeviceVector<int> iconserv;

    godunov::scheme godunov_scheme = godunov::scheme::PPM;
    std::string godunov_type;
    bool godunov_use_forces_in_trans{false};
    std::string advection_type{"Godunov"};
};

} // namespace amr_wind::pde

#endif /* ADVOP_GODUNOV_H */
