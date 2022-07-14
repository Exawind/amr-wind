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

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                for (int d = 0; d < ncomp; d++) {
                    farr(i, j, k, d) = x + y + z;
                }
            });
        }
    }
}

inline void init_field_mapped(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();
    const auto& nc_cc_fld = fld.repo().get_field("non_uniform_coord_cc");

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);
            const auto& cc_arr = nc_cc_fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = cc_arr(i, j, k, 0);
                const amrex::Real y = cc_arr(i, j, k, 1);
                const amrex::Real z = cc_arr(i, j, k, 2);

                for (int d = 0; d < ncomp; d++) {
                    farr(i, j, k, d) = x + y + z;
                }
            });
        }
    }
}

} // namespace amr_wind_tests

#endif /* TEST_ACT_UTILS_H */
