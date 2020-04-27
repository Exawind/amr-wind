#ifndef ICNS_OPS_H
#define ICNS_OPS_H

#include "PDETraits.H"
#include "PDEOps.H"
#include "PDEHelpers.H"
#include "AdvOp_Godunov.H"
#include "AdvOp_MOL.H"
#include "DiffusionOps.H"
#include "icns/icns.H"
#include "incflo_MAC_bcs.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace pde {

void advection_mac_project(FieldRepo& repo, FieldState fstate);

template<typename Scheme>
struct FieldRegOp<ICNS, Scheme>
{
    FieldRegOp(FieldRepo& repo_in) : repo(repo_in) {}

    PDEFields operator()(const SimTime& time, const int probtype)
    {
        auto fields = create_fields_instance<ICNS, Scheme>(time, repo, probtype);

        auto& rho = repo.declare_cc_field(
            "density", 1, Scheme::nghost_state, Scheme::num_states);
        auto& grad_p = repo.declare_cc_field("gp", ICNS::ndim, 0, 1);
        auto& pressure = repo.declare_nd_field("p", 1, 0, 1);
        repo.declare_face_normal_field(
            {"u_mac", "v_mac", "w_mac"}, 1, Scheme::nghost_mac, 1);

        rho.template register_fill_patch_op<FieldFillPatchOps<FieldBCDirichlet>>(
            repo.mesh(), time, probtype);
        grad_p.template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            repo.mesh(), time, probtype);
        pressure.template register_fill_patch_op<FieldFillConstScalar>(0.0);

        rho.fillpatch_on_regrid() = true;
        grad_p.fillpatch_on_regrid() = true;

        // Density always requires half state
        rho.create_state(FieldState::NPH);

        return fields;
    }

    FieldRepo& repo;
};

template<>
struct SrcTermOp<ICNS> : SrcTermOpBase<ICNS>
{
    SrcTermOp(PDEFields& fields_in)
    : SrcTermOpBase<ICNS>(fields_in)
    , grad_p(fields_in.repo.get_field("gp")) {}

    void operator()(const FieldState fstate)
    {
        const auto rhostate = field_impl::phi_state(fstate);
        auto& density = m_density.state(rhostate);

        const int nlevels = this->fields.repo.num_active_levels();
        for (int lev=0; lev < nlevels; ++lev) {
            auto& src_term = this->fields.src_term(lev);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(src_term, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& vf = src_term.array(mfi);
                const auto& rho = density(lev).const_array(mfi);
                const auto& gp = grad_p(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        amrex::Real rhoinv = 1.0 / rho(i, j, k);
                        vf(i, j, k, 0) = -(gp(i, j, k, 0)) * rhoinv;
                        vf(i, j, k, 1) = -(gp(i, j, k, 1)) * rhoinv;
                        vf(i, j, k, 2) = -(gp(i, j, k, 2)) * rhoinv;
                    });

                for (auto& src : this->sources) {
                    (*src)(lev, mfi, bx, fstate, vf);
                }
            }
        }
    }

    Field& grad_p;

};

template<>
struct AdvectionOp<ICNS, fvm::Godunov>
{
    AdvectionOp(PDEFields& fields_in)
    : fields(fields_in)
    , u_mac(fields_in.repo.get_field("u_mac"))
    , v_mac(fields_in.repo.get_field("v_mac"))
    , w_mac(fields_in.repo.get_field("w_mac"))
    {

        { // Prefix incflo
            amrex::ParmParse pp("incflo");
            pp.query("use_ppm", godunov_ppm);
            pp.query("godunov_use_forces_in_trans", godunov_use_forces_in_trans);
        }

        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(ICNS::ndim, 0);
    }

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();

        auto& src_term = fields.src_term;
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
            amrex::Box const& domain = geom[lev].Domain();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                amrex::FArrayBox scratch;
                for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    amrex::Box const& bx = mfi.tilebox();
                    amrex::Box const& bxg1 = amrex::grow(bx,1);
                    amrex::Box const& xbx = mfi.nodaltilebox(0);
                    amrex::Box const& ybx = mfi.nodaltilebox(1);
                    amrex::Box const& zbx = mfi.nodaltilebox(2);

                    amrex::Array4<amrex::Real> const& a_umac = u_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& a_vmac = v_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& a_wmac = w_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real const> const& a_vel = dof_field(lev).const_array(mfi);
                    amrex::Array4<amrex::Real const> const& a_f = src_term(lev).const_array(mfi);

                    scratch.resize(bxg1, ICNS::ndim*12+3);
    //                Elixir eli = scratch.elixir(); // not needed because of streamSynchronize later
                    amrex::Real* p = scratch.dataPtr();

                    amrex::Array4<amrex::Real> Imx = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Imx.size();
                    amrex::Array4<amrex::Real> Ipx = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Ipx.size();
                    amrex::Array4<amrex::Real> Imy = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Imy.size();
                    amrex::Array4<amrex::Real> Ipy = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Ipy.size();
                    amrex::Array4<amrex::Real> Imz = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Imz.size();
                    amrex::Array4<amrex::Real> Ipz = makeArray4(p,bxg1,ICNS::ndim);
                    p +=         Ipz.size();
                    amrex::Array4<amrex::Real> u_ad = makeArray4(p,amrex::Box(bx).grow(1,1).grow(2,1).surroundingNodes(0),1);
                    p +=         u_ad.size();
                    amrex::Array4<amrex::Real> v_ad = makeArray4(p,amrex::Box(bx).grow(0,1).grow(2,1).surroundingNodes(1),1);
                    p +=         v_ad.size();
                    amrex::Array4<amrex::Real> w_ad = makeArray4(p,amrex::Box(bx).grow(0,1).grow(1,1).surroundingNodes(2),1);
                    p +=         w_ad.size();

                    if (godunov_ppm){
                        godunov::predict_ppm (lev, bxg1, ICNS::ndim, Imx, Ipx, Imy, Ipy, Imz, Ipz, a_vel, a_vel,
                                              geom, dt, bcrec_device);
                    } else {
                        godunov::predict_plm_x (lev, bx, ICNS::ndim, Imx, Ipx, a_vel, a_vel,
                                     geom, dt, dof_field.bcrec(), bcrec_device);

                        godunov::predict_plm_y (lev, bx, ICNS::ndim, Imy, Ipy, a_vel, a_vel,
                                     geom, dt, dof_field.bcrec(), bcrec_device);

                        godunov::predict_plm_z (lev, bx, ICNS::ndim, Imz, Ipz, a_vel, a_vel,
                                     geom, dt, dof_field.bcrec(), bcrec_device);
                    }

                    godunov::make_trans_velocities(lev, amrex::Box(u_ad), amrex::Box(v_ad), amrex::Box(w_ad),
                                                   u_ad, v_ad, w_ad,
                                                   Imx, Ipx, Imy, Ipy, Imz, Ipz, a_vel, a_f,
                                                   geom, dt, bcrec_device, godunov_use_forces_in_trans);

                    godunov::predict_godunov(lev, bx, ICNS::ndim, xbx, ybx, zbx, a_umac, a_vmac, a_wmac,
                                             a_vel, u_ad, v_ad, w_ad, Imx, Ipx, Imy, Ipy, Imz, Ipz, a_f, p,
                                             geom, dt, bcrec_device, godunov_use_forces_in_trans);

                    incflo_set_mac_bcs(domain, xbx, ybx, zbx, a_umac, a_vmac, a_wmac, a_vel, dof_field.bcrec());

                    amrex::Gpu::streamSynchronize();  // otherwise we might be using too much memory
                }
            }
        }

        // MAC projection
        advection_mac_project(repo, fstate);

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
                amrex::FArrayBox tmpfab(
                    amrex::grow(bx, 1), ICNS::ndim * 14 + 1);

                godunov::compute_advection(
                    lev, bx, ICNS::ndim, conv_term(lev).array(mfi),
                    dof_field(lev).const_array(mfi),
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi), src_term(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), iconserv.data(),
                    tmpfab.dataPtr(), geom, dt, godunov_ppm);

                amrex::Gpu::streamSynchronize();
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    amrex::Gpu::DeviceVector<int> iconserv;

    bool godunov_ppm{true};
    bool godunov_use_forces_in_trans{false};
};

template<>
struct AdvectionOp<ICNS, fvm::MOL>
{
    AdvectionOp(PDEFields& fields_in)
    : fields(fields_in)
    , u_mac(fields_in.repo.get_field("u_mac"))
    , v_mac(fields_in.repo.get_field("v_mac"))
    , w_mac(fields_in.repo.get_field("w_mac"))
    {}

    void operator()(const FieldState fstate, const amrex::Real )
    {

        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();
        auto& conv_term = fields.conv_term.state(fstate);
        auto& dof_field = fields.field.state(fstate);

        //
        // Predict velocities
        //

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::Box const& domain = repo.mesh().Geom(lev).Domain();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            {
                for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    amrex::Box const& ubx = mfi.nodaltilebox(0);
                    amrex::Box const& vbx = mfi.nodaltilebox(1);
                    amrex::Box const& wbx = mfi.nodaltilebox(2);
                    amrex::Array4<amrex::Real> const& u = u_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& v = v_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real> const& w = w_mac(lev).array(mfi);
                    amrex::Array4<amrex::Real const> const& vcc = dof_field(lev).const_array(mfi);

                    mol::predict_vels_on_faces(lev,
                                               ubx, vbx, wbx,
                                               u, v, w,
                                               vcc,
                                               dof_field.bcrec(),
                                               dof_field.bcrec_device().data(),
                                               repo.mesh().Geom());

                    incflo_set_mac_bcs(domain,
                                       ubx, vbx, wbx,
                                       u, v, w,
                                       vcc,
                                       dof_field.bcrec());
                }
            }
        }

        advection_mac_project(repo, fstate);

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
};

template<>
struct TurbulenceOp<ICNS>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel,
                 PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        m_tmodel.update_mueff(m_fields.mueff);
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

}
}

#endif /* ICNS_OPS_H */
