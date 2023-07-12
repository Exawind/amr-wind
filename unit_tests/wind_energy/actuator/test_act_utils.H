#ifndef TEST_ACT_UTILS_H
#define TEST_ACT_UTILS_H

#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind_tests {

// Utility function to populate the velocity field used for tests
inline void init_field(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offsetx = 0.0;
    amrex::Real offsety = 0.0;
    amrex::Real offsetz = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offsetx = 0.5;
        offsety = 0.5;
        offsetz = 0.5;
    }
    if (fld.field_location() == amr_wind::FieldLoc::XFACE) {
        offsety = 0.5;
        offsetz = 0.5;
    }
    if (fld.field_location() == amr_wind::FieldLoc::YFACE) {
        offsetx = 0.5;
        offsetz = 0.5;
    }
    if (fld.field_location() == amr_wind::FieldLoc::ZFACE) {
        offsetx = 0.5;
        offsety = 0.5;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + offsetx) * dx[0];
                const amrex::Real y = problo[1] + (j + offsety) * dx[1];
                const amrex::Real z = problo[2] + (k + offsetz) * dx[2];

                for (int d = 0; d < ncomp; d++) {
                    farr(i, j, k, d) = x + y + z;
                }
            });
        }
    }
}

} // namespace amr_wind_tests

#endif /* TEST_ACT_UTILS_H */
