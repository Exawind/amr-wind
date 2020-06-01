#ifndef MLMGOPTIONS_H
#define MLMGOPTIONS_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

/** Collection of solver options to control MLMG behavior
 */
struct MLMGOptions
{
    /** Parse user options with a given prefix
     *
     *  \param prefix Prefix used to parse user inputs, e.g., mac_proj
     */
    MLMGOptions(const std::string& prefix) { parse_options(prefix); }

    /** Parse options in two stages
     *
     *  This constructor parses the default options from `default_prefix` and
     *  then overrides certain options based on a custom prefix.
     *
     *  \param default_prefix Namespace for default options
     *  \param custom_prefix Namespace for override options
     */
    MLMGOptions(
        const std::string& default_prefix, const std::string& custom_prefix)
    {
        parse_options(default_prefix);
        parse_options(custom_prefix);
    }

    //! Bottom solver type bicgstab, cg, hypre, etc.
    std::string bottom_solver_type{"bicgstab"};

    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real rel_tol{1.0e-11};
    //! Absolute tolerance for convergence checks
    amrex::Real abs_tol{1.0e-14};

    int verbose{0};
    int max_iter{200};
    int max_coarsen_level{100};
    int max_order{2};
    int fmg_max_iter{0};

    int cg_verbose{0};
    int cg_max_iter{200};

private:
    void parse_options(const std::string& prefix)
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
};

} // namespace amr_wind

#endif /* MLMGOPTIONS_H */
