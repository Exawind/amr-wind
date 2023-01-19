#include "amr-wind/core/IntScratchField.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Gpu.H"
#include "AMReX_IArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_iMultiFab.H"

namespace amr_wind {


namespace {


} // namespace


void IntScratchField::setVal(int value) noexcept
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 1");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value);
    }
}

/*
void IntScratchField::setVal(
    int value, int start_comp, int num_comp, int nghost) noexcept
{
    BL_PROFILE("amr-wind::IntScratchField::setVal 2");
    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        operator()(lev).setVal(value, start_comp, num_comp, nghost);
    }
}
*/
/*
void IntScratchField::setVal(const amrex::Vector<int>& values, int nghost) noexcept
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
*/


} // namespace amr_wind
