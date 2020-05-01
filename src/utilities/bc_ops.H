#ifndef BC_OPS_H
#define BC_OPS_H

/** \file bc_ops.H
 *
 *  Boundary condition utilities
 */

#include <AMReX_BCRec.H>

namespace amr_wind {
namespace utils {

//! Return a pair of bools,
//! value is true if at least one boundary has an ext_dir specified
std::pair<bool, bool> has_extdir(amrex::BCRec const* bcrec, int ncomp, int dir);

} // namespace utils
} // namespace amr_wind

#endif /* BC_OPS_H */
