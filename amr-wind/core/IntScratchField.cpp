#include "amr-wind/core/IntScratchField.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Gpu.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_iMultiFab.H"

namespace amr_wind {


namespace {
struct ISFBCNoOp
{
    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& /* iv */,
        amrex::Array4<int> const& /* field */,
        const int /* dcomp */,
        const int /* numcomp */,
        amrex::GeometryData const& /* geom */,
        const amrex::Real /* time */,
        const amrex::BCRec* /* bcr */,
        const int /* bcomp */,
        const int /* orig_comp */) const
    {}
};

amrex::Vector<amrex::BCRec>
scratch_field_bcrec(const amrex::Geometry& geom, const int ncomp)
{
    amrex::Vector<amrex::BCRec> bcrec(ncomp);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto bcrv = geom.isPeriodic(dir) ? amrex::BCType::int_dir
                                         : amrex::BCType::hoextrap;

        for (int n = 0; n < ncomp; ++n) {
            bcrec[n].setLo(dir, bcrv);
            bcrec[n].setHi(dir, bcrv);
        }
    }

    return bcrec;
}

} // namespace
void IntScratchField::fillpatch(amrex::Real time) noexcept
{
    fillpatch(time, num_grow());
}

void IntScratchField::fillpatch(
    amrex::Real time, const amrex::IntVect& ng) noexcept
{
    const int nlevels = repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        fillpatch(lev, time, this->operator()(lev), ng);
    }
}

void IntScratchField::fillpatch(
    int lev,
    amrex::Real time,
    amrex::iMultiFab& mfab,
    const amrex::IntVect& nghost) noexcept
{
    const auto& mesh = repo().mesh();
    auto bcrec = scratch_field_bcrec(mesh.Geom(lev), num_comp());
    fillpatch(lev, time, mfab, nghost, bcrec);
}

void IntScratchField::fillpatch(
    int lev,
    amrex::Real time,
    amrex::iMultiFab& mfab,
    const amrex::IntVect& nghost,
    amrex::Vector<amrex::BCRec>& bcrec) noexcept
{
    const auto& mesh = repo().mesh();
    amrex::Interpolater* mapper = &amrex::cell_cons_interp;

    if (lev == 0) {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<ISFBCNoOp>> physbc(
            mesh.Geom(lev), bcrec, ISFBCNoOp{});

        amrex::FillPatchSingleLevel(
            mfab, nghost, time, {&this->operator()(lev)}, {time}, 0, 0,
            num_comp(), mesh.Geom(lev), physbc, 0);
    } else {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<ISFBCNoOp>> cphysbc(
            mesh.Geom(lev - 1), bcrec, ISFBCNoOp{});

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<ISFBCNoOp>> fphysbc(
            mesh.Geom(lev), bcrec, ISFBCNoOp{});

        amrex::FillPatchTwoLevels(
            mfab, nghost, time, {&this->operator()(lev - 1)}, {time},
            {&this->operator()(lev)}, {time}, 0, 0, num_comp(),
            mesh.Geom(lev - 1), mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
            mesh.refRatio(lev - 1), mapper, bcrec, 0);
    }
}

} // namespace amr_wind
