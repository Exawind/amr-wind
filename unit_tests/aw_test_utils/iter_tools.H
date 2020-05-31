#ifndef ITER_TOOLS_H
#define ITER_TOOLS_H

#include "AMReX_MultiFab.H"
#include "AMReX_Vector.H"
#include "amr-wind/core/FieldRepo.H"

/** \file iter_tools.H
 *  Iteration utilities
 */

namespace amr_wind_tests {

/** Helper to run a functor over all levels of a given field
 *
 *  The functor is a void function that accepts two arguments: the current
 *  level, and an MFIter instance for the MultiFab at the given level.
 *
 *  @param nlevels Number of valid levels in the AMR hierarchy
 *  @param field MultiFabs for each of the `nlevels` levels
 *  @param func Functor to execute over the MultiFabs
 */
template <typename Functor>
void run_algorithm(
    const int nlevels, amrex::Vector<amrex::MultiFab*>& field, Functor func)
{
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& lfab = *field[lev];

        for (amrex::MFIter mfi(lfab); mfi.isValid(); ++mfi) {
            func(lev, mfi);
        }
    }
}

template <typename Functor>
void run_algorithm(amr_wind::Field& field, Functor func)
{
    const int nlevels = field.repo().num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        auto& lfab = field(lev);

        for (amrex::MFIter mfi(lfab); mfi.isValid(); ++mfi) {
            func(lev, mfi);
        }
    }
}

} // namespace amr_wind_tests

#endif /* ITER_TOOLS_H */
