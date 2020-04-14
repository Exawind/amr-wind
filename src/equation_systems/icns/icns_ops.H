#ifndef ICNS_OPS_H
#define ICNS_OPS_H

#include "PDETraits.H"
#include "PDEOps.H"
#include "PDEHelpers.H"
#include "icns/icns.H"

namespace amr_wind {
namespace pde {

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

        return fields;
    }

    FieldRepo& repo;
};

template<>
struct SrcTermOp<ICNS> : SrcTermOpBase<ICNS>
{
    SrcTermOp(PDEFields& fields_in) : SrcTermOpBase<ICNS>(fields_in) {}

    void operator()(const FieldState fstate)
    {
        const auto rhostate = field_impl::phi_state(fstate);
        auto& density = fields.repo.get_field("density", rhostate);
        auto& grad_p = fields.repo.get_field("gp");

        const int nlevels = this->fields.repo.num_active_levels();
        for (int lev=0; lev < nlevels; ++lev) {
            auto& src_term = this->fields.src_term(lev);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
};

}
}

#endif /* ICNS_OPS_H */
