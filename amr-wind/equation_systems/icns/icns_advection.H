#ifndef ICNS_ADVECTION_H
#define ICNS_ADVECTION_H

#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_MacProjector.H"

namespace amr_wind {
namespace pde {

class MacProjOp
{
public:
    using FaceFabPtrVec =
        amrex::Vector<amrex::Array<const amrex::MultiFab*, ICNS::ndim>>;

    MacProjOp(FieldRepo&, bool, bool);

    void operator()(const FieldState fstate, const amrex::Real dt);

private:
    void init_projector(const FaceFabPtrVec&) noexcept;
    void init_projector(const amrex::Real) noexcept;

    FieldRepo& m_repo;
    std::unique_ptr<amrex::MacProjector> m_mac_proj;
    MLMGOptions m_options;
    bool m_has_overset{false};
    bool m_need_init{true};
    bool m_variable_density{false};
    amrex::Real m_rho_0{1.0};
};

/** Godunov scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::Godunov>
{
    AdvectionOp(PDEFields& fields_in, bool has_overset, bool variable_density)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_macproj_op(fields.repo, has_overset, variable_density)
    {

        amrex::ParmParse pp("incflo");
        pp.query("godunov_type", godunov_type);
        if (pp.contains("use_ppm") || pp.contains("use_limiter"))
            amrex::Abort(
                "Godunov: use_ppm and use_limiter are deprecated. Please "
                "update input file");

        if (amrex::toLower(godunov_type) == "plm") {
            godunov_scheme = godunov::scheme::PLM;
        } else if (amrex::toLower(godunov_type) == "ppm") {
            godunov_scheme = godunov::scheme::PPM;
        } else if (amrex::toLower(godunov_type) == "ppm_nolim") {
            godunov_scheme = godunov::scheme::PPM_NOLIM;
        } else if (amrex::toLower(godunov_type) == "weno") {
            godunov_scheme = godunov::scheme::WENO;
        } else {
            amrex::Print() << "For godunov_type select between plm, ppm, "
                              "ppm_nolim and weno: it defaults to ppm"
                           << std::endl;
            godunov_scheme = godunov::scheme::PPM;
        }

        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(ICNS::ndim, 0);
    }

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        auto& conv_term = fields.conv_term;
        auto& dof_field = fields.field.state(fstate);

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
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                amrex::FArrayBox scratch;
                for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU());
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
                    //                needed because of streamSynchronize later
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
                        amrex::Box(bx).grow(1, 1).grow(2, 1).surroundingNodes(
                            0),
                        1);
                    p += u_ad.size();
                    amrex::Array4<amrex::Real> v_ad = makeArray4(
                        p,
                        amrex::Box(bx).grow(0, 1).grow(2, 1).surroundingNodes(
                            1),
                        1);
                    p += v_ad.size();
                    amrex::Array4<amrex::Real> w_ad = makeArray4(
                        p,
                        amrex::Box(bx).grow(0, 1).grow(1, 1).surroundingNodes(
                            2),
                        1);
                    p += w_ad.size();

                    switch (godunov_scheme) {
                    case godunov::scheme::PPM: {
                        godunov::predict_ppm(
                            lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz, Ipz,
                            a_vel, a_vel, geom, dt, bcrec_device, true);
                        break;
                    }
                    case godunov::scheme::PPM_NOLIM: {
                        godunov::predict_ppm(
                            lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz, Ipz,
                            a_vel, a_vel, geom, dt, bcrec_device, false);
                        break;
                    }
                    case godunov::scheme::WENO: {
                        godunov::predict_weno(
                            lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz, Ipz,
                            a_vel, a_vel, geom, dt, bcrec_device);
                        break;
                    }
                    case godunov::scheme::PLM: {
                        godunov::predict_plm_x(
                            lev, bx, ICNS::ndim, Imx, Ipx, a_vel, a_vel, geom,
                            dt, dof_field.bcrec(), bcrec_device);

                        godunov::predict_plm_y(
                            lev, bx, ICNS::ndim, Imy, Ipy, a_vel, a_vel, geom,
                            dt, dof_field.bcrec(), bcrec_device);

                        godunov::predict_plm_z(
                            lev, bx, ICNS::ndim, Imz, Ipz, a_vel, a_vel, geom,
                            dt, dof_field.bcrec(), bcrec_device);
                        break;
                    }
                    }

                    godunov::make_trans_velocities(
                        lev, amrex::Box(u_ad), amrex::Box(v_ad),
                        amrex::Box(w_ad), u_ad, v_ad, w_ad, Imx, Ipx, Imy, Ipy,
                        Imz, Ipz, a_vel, a_f, geom, dt, bcrec_device,
                        godunov_use_forces_in_trans);

                    godunov::predict_godunov(
                        lev, bx, ICNS::ndim, xbx, ybx, zbx, a_umac, a_vmac,
                        a_wmac, a_vel, u_ad, v_ad, w_ad, Imx, Ipx, Imy, Ipy,
                        Imz, Ipz, a_f, p, geom, dt, bcrec_device,
                        godunov_use_forces_in_trans);

                    amrex::Gpu::streamSynchronize(); // otherwise we might be
                                                     // using too much memory
                }
            }
        }

        // MAC projection
        m_macproj_op(fstate, dt);

        //
        // Advect momentum eqns
        //
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            u_mac(lev).FillBoundary(geom[lev].periodicity());
            v_mac(lev).FillBoundary(geom[lev].periodicity());
            w_mac(lev).FillBoundary(geom[lev].periodicity());

            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion())
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), ICNS::ndim * 14);

                godunov::compute_advection(
                    lev, bx, ICNS::ndim, conv_term(lev).array(mfi),
                    dof_field(lev).const_array(mfi),
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi), src_term(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), iconserv.data(),
                    tmpfab.dataPtr(), geom, dt, godunov_scheme);

                amrex::Gpu::streamSynchronize();
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
    std::string godunov_type;
    bool godunov_use_forces_in_trans{false};
};

/** MOL scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::MOL>
{
    AdvectionOp(PDEFields& fields_in, bool has_overset, bool variable_density)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_macproj_op(fields.repo, has_overset, variable_density)
    {}

    void operator()(const FieldState fstate, const amrex::Real dt)
    {

        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();
        auto& conv_term = fields.conv_term.state(fstate);
        auto& dof_field = fields.field.state(fstate);

        //
        // Predict velocities
        //

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi) {
                    amrex::Box const& ubx = mfi.nodaltilebox(0);
                    amrex::Box const& vbx = mfi.nodaltilebox(1);
                    amrex::Box const& wbx = mfi.nodaltilebox(2);
                    amrex::Array4<amrex::Real> const& u = u_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& v = v_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& w = w_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real const> const& vcc =
                        dof_field(lev).const_array(mfi);

                    mol::predict_vels_on_faces(
                        lev, ubx, vbx, wbx, u, v, w, vcc, dof_field.bcrec(),
                        dof_field.bcrec_device().data(), repo.mesh().Geom());
                }
            }
        }

        m_macproj_op(fstate, dt);

        //
        // Advect velocity
        //

        int nmaxcomp = AMREX_SPACEDIM;
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            // if (amrex::Gpu::notInLaunchRegion())
            // mfi_info.EnableTiling(amrex::IntVect(1024,16,16)).SetDynamic(true);
            if (amrex::Gpu::notInLaunchRegion())
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                amrex::Box const& bx = mfi.tilebox();

                amrex::Box tmpbox = amrex::surroundingNodes(bx);
                const int tmpcomp = nmaxcomp * AMREX_SPACEDIM;

                amrex::FArrayBox tmpfab(tmpbox, tmpcomp);
                amrex::Elixir eli = tmpfab.elixir();

                amrex::Array4<amrex::Real> fx = tmpfab.array(nmaxcomp * 0);
                amrex::Array4<amrex::Real> fy = tmpfab.array(nmaxcomp * 1);
                amrex::Array4<amrex::Real> fz = tmpfab.array(nmaxcomp * 2);

                mol::compute_convective_fluxes(
                    lev, bx, AMREX_SPACEDIM, fx, fy, fz,
                    dof_field(lev).const_array(mfi),
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

    MacProjOp m_macproj_op;
};

} // namespace pde
} // namespace amr_wind

#endif /* ICNS_ADVECTION_H */
