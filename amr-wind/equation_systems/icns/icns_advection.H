#ifndef ICNS_ADVECTION_H
#define ICNS_ADVECTION_H

#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/vof/vof_momentum_flux.H"

#include "AMReX_MultiFabUtil.H"
#include "hydro_MacProjector.H"
#include "hydro_mol.H"
#include "hydro_utils.H"

namespace amr_wind::pde {

class MacProjOp
{
public:
    using FaceFabPtrVec =
        amrex::Vector<amrex::Array<const amrex::MultiFab*, ICNS::ndim>>;

    MacProjOp(
        FieldRepo& /*repo*/,
        bool /*has_overset*/,
        bool /*variable_density*/,
        bool /*mesh_mapping*/);

    void operator()(const FieldState fstate, const amrex::Real dt);

    static void mac_proj_to_uniform_space(
        const amr_wind::FieldRepo& /*repo*/,
        amr_wind::Field& /*u_mac*/,
        amr_wind::Field& /*v_mac*/,
        amr_wind::Field& /*w_mac*/,
        amrex::Array<amrex::MultiFab*, ICNS::ndim>& /*rho_face*/,
        amrex::Real /*ovst_fac*/,
        int /*lev*/) noexcept;

    amrex::Real rho0() const { return m_rho_0; }

private:
    void init_projector(const FaceFabPtrVec& /*beta*/) noexcept;
    void init_projector(const amrex::Real /*beta*/) noexcept;

    FieldRepo& m_repo;
    std::unique_ptr<Hydro::MacProjector> m_mac_proj;
    MLMGOptions m_options;
    bool m_has_overset{false};
    bool m_need_init{true};
    bool m_variable_density{false};
    bool m_mesh_mapping{false};
    amrex::Real m_rho_0{1.0};
};

/** Godunov scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::Godunov>
{
    AdvectionOp(
        PDEFields& fields_in,
        bool has_overset,
        bool variable_density,
        bool mesh_mapping)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_macproj_op(fields.repo, has_overset, variable_density, mesh_mapping)
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
            // use Godunov for premac, use BDS for postmac. Eventually
            // there will be a premac BDS
            postmac_advection_type = "BDS";
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

        // Flux calculation used in multiphase portions of domain
        pp.query("mflux_type", mflux_type);
        if (amrex::toLower(mflux_type) == "minmod") {
            mflux_scheme = godunov::scheme::MINMOD;
        } else if (amrex::toLower(mflux_type) == "upwind") {
            mflux_scheme = godunov::scheme::UPWIND;
        } else {
            amrex::Print() << "For mflux_type, default is upwind" << std::endl;
            mflux_scheme = godunov::scheme::UPWIND;
        }

        // Formulation of discrete ICNS equation
        // 1 = conservative (default), 0 = nonconservative
        pp.query("icns_conserv", m_cons);
        iconserv.resize(ICNS::ndim, m_cons);
    }

    void preadvect(const FieldState fstate, const amrex::Real dt)
    {

        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        const auto& dof_field = fields.field.state(fstate);
        auto bcrec_device = dof_field.bcrec_device();

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            u_mac(lev).setBndry(0.0);
            v_mac(lev).setBndry(0.0);
            w_mac(lev).setBndry(0.0);
        }

        //
        // Predict
        //
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            if ((godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                (godunov_scheme == godunov::scheme::WENOJS) ||
                (godunov_scheme == godunov::scheme::WENOZ)) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                {
                    amrex::FArrayBox scratch;
                    for (amrex::MFIter mfi(
                             dof_field(lev), amrex::TilingIfNotGPU());
                         mfi.isValid(); ++mfi) {
                        amrex::Box const& bx = mfi.tilebox();
                        amrex::Box const& bxg1 = amrex::grow(bx, 1);
                        amrex::Box const& xbx = mfi.nodaltilebox(0);
                        amrex::Box const& ybx = mfi.nodaltilebox(1);
                        amrex::Box const& zbx = mfi.nodaltilebox(2);

                        amrex::Array4<amrex::Real> const& a_umac =
                            u_mac(lev).array(mfi);
                        amrex::Array4<amrex::Real> const& a_vmac =
                            v_mac(lev).array(mfi);
                        amrex::Array4<amrex::Real> const& a_wmac =
                            w_mac(lev).array(mfi);
                        amrex::Array4<amrex::Real const> const& a_vel =
                            dof_field(lev).const_array(mfi);
                        amrex::Array4<amrex::Real const> const& a_f =
                            src_term(lev).const_array(mfi);

                        scratch.resize(bxg1, ICNS::ndim * 12 + 3);
                        //                Elixir eli = scratch.elixir(); // not
                        //                needed because of streamSynchronize
                        //                later
                        amrex::Real* p = scratch.dataPtr();

                        amrex::Array4<amrex::Real> Imx =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Imx.size();
                        amrex::Array4<amrex::Real> Ipx =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Ipx.size();
                        amrex::Array4<amrex::Real> Imy =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Imy.size();
                        amrex::Array4<amrex::Real> Ipy =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Ipy.size();
                        amrex::Array4<amrex::Real> Imz =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Imz.size();
                        amrex::Array4<amrex::Real> Ipz =
                            makeArray4(p, bxg1, ICNS::ndim);
                        p += Ipz.size();
                        amrex::Array4<amrex::Real> u_ad = makeArray4(
                            p,
                            amrex::Box(bx)
                                .grow(1, 1)
                                .grow(2, 1)
                                .surroundingNodes(0),
                            1);
                        p += u_ad.size();
                        amrex::Array4<amrex::Real> v_ad = makeArray4(
                            p,
                            amrex::Box(bx)
                                .grow(0, 1)
                                .grow(2, 1)
                                .surroundingNodes(1),
                            1);
                        p += v_ad.size();
                        amrex::Array4<amrex::Real> w_ad = makeArray4(
                            p,
                            amrex::Box(bx)
                                .grow(0, 1)
                                .grow(1, 1)
                                .surroundingNodes(2),
                            1);
                        p += w_ad.size();

                        switch (godunov_scheme) {
                        case godunov::scheme::PPM_NOLIM: {
                            godunov::predict_ppm(
                                lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz,
                                Ipz, a_vel, a_vel, geom, dt, bcrec_device,
                                false);
                            break;
                        }
                        case godunov::scheme::WENOJS: {
                            godunov::predict_weno(
                                lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz,
                                Ipz, a_vel, a_vel, geom, dt, bcrec_device,
                                true);
                            break;
                        }
                        case godunov::scheme::WENOZ: {
                            godunov::predict_weno(
                                lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz,
                                Ipz, a_vel, a_vel, geom, dt, bcrec_device,
                                false);
                            break;
                        }
                        default: {
                            amrex::Abort(
                                "Only PPM_NOLIM, WENOZ, and WENOJS use this "
                                "code path");
                        }
                        }

                        godunov::make_trans_velocities(
                            lev, amrex::Box(u_ad), amrex::Box(v_ad),
                            amrex::Box(w_ad), u_ad, v_ad, w_ad, Imx, Ipx, Imy,
                            Ipy, Imz, Ipz, a_vel, a_f, geom, dt, bcrec_device,
                            godunov_use_forces_in_trans);

                        godunov::predict_godunov(
                            lev, bx, ICNS::ndim, xbx, ybx, zbx, a_umac, a_vmac,
                            a_wmac, a_vel, u_ad, v_ad, w_ad, Imx, Ipx, Imy, Ipy,
                            Imz, Ipz, a_f, p, geom, dt, bcrec_device,
                            godunov_use_forces_in_trans);

                        amrex::Gpu::streamSynchronize(); // otherwise we might
                                                         // be using too much
                                                         // memory
                    }
                }
            } else if (
                (godunov_scheme == godunov::scheme::PPM) ||
                (godunov_scheme == godunov::scheme::PLM) ||
                (godunov_scheme == godunov::scheme::BDS)) {
                const bool godunov_use_ppm =
                    godunov_scheme == godunov::scheme::PPM;
                HydroUtils::ExtrapVelToFaces(
                    dof_field(lev), src_term(lev), u_mac(lev), v_mac(lev),
                    w_mac(lev), dof_field.bcrec(),
                    dof_field.bcrec_device().data(), repo.mesh().Geom(lev), dt,
                    godunov_use_ppm, godunov_use_forces_in_trans,
                    premac_advection_type);
            } else {
                amrex::Abort("Invalid godunov scheme");
            }
        }

        // MAC projection
        m_macproj_op(fstate, dt);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            u_mac(lev).FillBoundary(geom[lev].periodicity());
            v_mac(lev).FillBoundary(geom[lev].periodicity());
            w_mac(lev).FillBoundary(geom[lev].periodicity());
        }
    }

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        // cppcheck-suppress constVariable
        auto& conv_term = fields.conv_term;
        const auto& dof_field = fields.field.state(fstate);

        auto flux_x =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::ZFACE);
        auto face_x =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto face_y =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto face_z =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::ZFACE);

        const auto& rho_o =
            repo.get_field("density").state(amr_wind::FieldState::Old);

        const bool mphase_vof = repo.field_exists("vof");

        //
        // Advect momentum eqns
        //
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {

            // form multifab for transport variable and source term
            amrex::MultiFab q(
                dof_field(lev).boxArray(), dof_field(lev).DistributionMap(),
                ICNS::ndim, fvm::Godunov::nghost_state);
            amrex::MultiFab::Copy(
                q, dof_field(lev), 0, 0, ICNS::ndim,
                fvm::Godunov::nghost_state);
            amrex::MultiFab fq(
                src_term(lev).boxArray(), src_term(lev).DistributionMap(),
                ICNS::ndim, fvm::Godunov::nghost_src);
            amrex::MultiFab::Copy(
                fq, src_term(lev), 0, 0, ICNS::ndim, fvm::Godunov::nghost_src);

            // Calculate fluxes using momentum directly
            if (!mphase_vof) {
                amrex::MultiFab::Multiply(
                    q, rho_o(lev), 0, 0, 1, fvm::Godunov::nghost_state);
                amrex::MultiFab::Multiply(
                    q, rho_o(lev), 0, 1, 1, fvm::Godunov::nghost_state);
                amrex::MultiFab::Multiply(
                    q, rho_o(lev), 0, 2, 1, fvm::Godunov::nghost_state);
                // Source terms are at old state during calculation of advection
                // terms
                amrex::MultiFab::Multiply(
                    fq, rho_o(lev), 0, 0, 1, fvm::Godunov::nghost_src);
                amrex::MultiFab::Multiply(
                    fq, rho_o(lev), 0, 1, 1, fvm::Godunov::nghost_src);
                amrex::MultiFab::Multiply(
                    fq, rho_o(lev), 0, 2, 1, fvm::Godunov::nghost_src);
            }

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

                // These schemes are not in amrex-hydro yet
                if ((godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                    (godunov_scheme == godunov::scheme::WENOJS) ||
                    (godunov_scheme == godunov::scheme::WENOZ)) {
                    amrex::FArrayBox tmpfab(
                        amrex::grow(bx, 1), ICNS::ndim * 14);
                    godunov::compute_fluxes(
                        lev, bx, ICNS::ndim, (*flux_x)(lev).array(mfi),
                        (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                        q.const_array(mfi), u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi), fq.const_array(mfi),
                        dof_field.bcrec_device().data(), iconserv.data(),
                        tmpfab.dataPtr(), geom, dt, godunov_scheme);
                } else if (
                    (godunov_scheme == godunov::scheme::PPM) ||
                    (godunov_scheme == godunov::scheme::PLM) ||
                    (godunov_scheme == godunov::scheme::BDS)) {
                    amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 1);
                    tmpfab.setVal<amrex::RunOn::Device>(0.0);
                    const auto& divu = tmpfab.array();
                    const bool is_velocity = true;
                    const bool known_edge_state = false;
                    const bool godunov_use_ppm =
                        godunov_scheme == godunov::scheme::PPM;
                    const bool fluxes_are_area_weighted = false;
                    HydroUtils::ComputeFluxesOnBoxFromState(
                        bx, ICNS::ndim, mfi, q.const_array(mfi),
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
                        divu, fq.const_array(mfi), geom[lev], dt,
                        dof_field.bcrec(), dof_field.bcrec_device().data(),
                        iconserv.data(), godunov_use_ppm,
                        godunov_use_forces_in_trans, is_velocity,
                        fluxes_are_area_weighted, postmac_advection_type);
                } else {
                    amrex::Abort("Invalid godunov scheme");
                }

                amrex::Gpu::streamSynchronize();
            }
        }

        // Multiphase flux operations
        if (mphase_vof) {
            // Loop levels
            multiphase::hybrid_fluxes(
                repo, ICNS::ndim, iconserv, (*flux_x), (*flux_y), (*flux_z),
                dof_field, src_term, rho_o, u_mac, v_mac, w_mac,
                dof_field.bcrec_device().data(), rho_o.bcrec_device().data(),
                dt, mflux_scheme);
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
                    lev, bx, ICNS::ndim, conv_term(lev).array(mfi),
                    (*flux_x)(lev).array(mfi), (*flux_y)(lev).array(mfi),
                    (*flux_z)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    iconserv.data(), geom);
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;

    MacProjOp m_macproj_op;
    amrex::Gpu::DeviceVector<int> iconserv;

    godunov::scheme godunov_scheme = godunov::scheme::PPM;
    godunov::scheme mflux_scheme = godunov::scheme::UPWIND;
    std::string godunov_type;
    std::string mflux_type;
    bool godunov_use_forces_in_trans{false};
    int m_cons{1};
    std::string premac_advection_type{"Godunov"};
    std::string postmac_advection_type{"Godunov"};
};

/** MOL scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::MOL>
{
    AdvectionOp(
        PDEFields& fields_in,
        bool has_overset,
        bool variable_density,
        bool mesh_mapping)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_mesh_mapping(mesh_mapping)
        , m_macproj_op(
              fields.repo, has_overset, variable_density, m_mesh_mapping)
    {}

    void preadvect(const FieldState fstate, const amrex::Real dt)
    {

        // cppcheck-suppress constVariable
        auto& repo = fields.repo;
        auto& dof_field = fields.field.state(fstate);

        // computation of velocity on faces requires
        // dof field to be in stretched mesh space
        if (dof_field.in_uniform_space() && m_mesh_mapping) {
            dof_field.to_stretched_space();
        }

        //
        // Predict velocities
        //

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            MOL::ExtrapVelToFaces(
                dof_field(lev), u_mac(lev), v_mac(lev), w_mac(lev),
                repo.mesh().Geom(lev), dof_field.bcrec(),
                dof_field.bcrec_device().data());
        }

        m_macproj_op(fstate, dt);
    }

    void operator()(const FieldState fstate, const amrex::Real /*unused*/)
    {

        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();
        // cppcheck-suppress constVariable
        auto& conv_term = fields.conv_term.state(fstate);
        const auto& dof_field = fields.field.state(fstate);
        const auto& rho = repo.get_field("density").state(fstate);

        //
        // Advect velocity
        //

        int nmaxcomp = AMREX_SPACEDIM;
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {

            amrex::MFItInfo mfi_info;
            // if (amrex::Gpu::notInLaunchRegion())
            // mfi_info.EnableTiling(amrex::IntVect(1024,16,16)).SetDynamic(true);
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                amrex::Box const& bx = mfi.tilebox();
                amrex::Box gbx = grow(bx, fvm::MOL::nghost_state);

                // Set up momentum array
                amrex::FArrayBox qfab(gbx, ICNS::ndim);
                const auto& q = qfab.array();
                // Calculate momentum
                auto rho_arr = rho(lev).const_array(mfi);
                auto vel_arr = dof_field(lev).const_array(mfi);
                amrex::ParallelFor(
                    gbx, ICNS::ndim,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        q(i, j, k, n) = rho_arr(i, j, k) * vel_arr(i, j, k, n);
                    });
                // Doing this explicitly, instead of through a Multiply command,
                // helps avoid floating-point errors with intel compilers and
                // mimics the implementation in equation_systems/AdvOp_MOL.H
                amrex::Gpu::streamSynchronize();

                amrex::Box tmpbox = amrex::surroundingNodes(bx);
                const int tmpcomp = nmaxcomp * AMREX_SPACEDIM;

                amrex::FArrayBox tmpfab(tmpbox, tmpcomp);
                amrex::Elixir eli = tmpfab.elixir();

                amrex::Array4<amrex::Real> fx = tmpfab.array(0);
                amrex::Array4<amrex::Real> fy = tmpfab.array(nmaxcomp);
                amrex::Array4<amrex::Real> fz = tmpfab.array(nmaxcomp * 2);

                mol::compute_convective_fluxes(
                    lev, bx, AMREX_SPACEDIM, fx, fy, fz, q,
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi), dof_field.bcrec().data(),
                    dof_field.bcrec_device().data(), geom);

                mol::compute_convective_rate(
                    bx, AMREX_SPACEDIM, conv_term(lev).array(mfi), fx, fy, fz,
                    geom[lev].InvCellSizeArray());
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;

    bool m_mesh_mapping;

    MacProjOp m_macproj_op;
};

} // namespace amr_wind::pde

#endif /* ICNS_ADVECTION_H */
