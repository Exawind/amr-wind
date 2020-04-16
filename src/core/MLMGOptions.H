#ifndef MLMGOPTIONS_H
#define MLMGOPTIONS_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

struct MLMGOptions
{
    MLMGOptions(const std::string& prefix)
    {
        amrex::ParmParse pp(prefix);

        pp.query("mg_verbose", verbose);
        pp.query("mg_cg_verbose", cg_verbose);

        pp.query("mg_max_iter", max_iter);
        pp.query("mg_cg_max_iter", cg_max_iter);
        pp.query("mg_fmg_max_iter", fmg_max_iter);
        pp.query("mg_max_coarsening_level", max_coarsen_level);
        pp.query("mg_max_order", max_order);
        pp.query("mg_rtol", rel_tol);
        pp.query("mg_atol", abs_tol);
        pp.query("bottom_solver_type", bottom_solver_type);
    }

    std::string bottom_solver_type{"bicgstab"};
    amrex::Real rel_tol{1.0e-11};
    amrex::Real abs_tol{1.0e-14};

    int verbose{0};
    int max_iter{200};
    int max_coarsen_level{100};
    int max_order{2};
    int fmg_max_iter{0};

    int cg_verbose{0};
    int cg_max_iter{200};
};

}


#endif /* MLMGOPTIONS_H */
