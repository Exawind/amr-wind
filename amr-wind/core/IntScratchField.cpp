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

/*
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
*/

} // namespace


void IntScratchField::setVal(int value) noexcept
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 1");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value);
    }
}

void IntScratchField::setVal(
    int value, int start_comp, int num_comp, int nghost) noexcept
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 2");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value, start_comp, num_comp, nghost);
    }
}

void IntField::setVal(const amrex::Vector<int>& values, int nghost) noexcept
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 3");
    AMREX_ASSERT(num_comp() == static_cast<int>(values.size()));

    // Update 1 component at a time
    const int ncomp = 1;
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        auto& mf = operator()(lev);
        for (int ic = 0; ic < num_comp(); ++ic) {
            int value = values[ic];
            mf.setVal(value, ic, ncomp, nghost);
        }
    }
}

/*
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
*/

} // namespace amr_wind
