#ifndef MLMGOPTIONS_H
#define MLMGOPTIONS_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MLLinOp.H"
#include "AMReX_MLMG.H"

namespace Hydro {
class NodalProjector;
class MacProjector;
} // namespace Hydro

namespace amr_wind {

/** Interface to control the behavior of AMReX LinearSolvers
 *
 *  MLMGOptions provides a unified interface to set options to the linear
 *  operator object, MLMG solver, external solvers (e.g., hypre), as well as
 *  MacProjector and NodalProjector objects.
 */
struct MLMGOptions
{
    /** Parse user options with a given prefix
     *
     *  \param prefix Prefix used to parse user inputs, e.g., mac_proj
     */
    explicit MLMGOptions(const std::string& prefix);

    /** Parse options in two stages
     *
     *  This constructor parses the default options from `default_prefix` and
     *  then overrides certain options based on a custom prefix.
     *
     *  \param default_prefix Namespace for default options
     *  \param custom_prefix Namespace for override options
     */
    MLMGOptions(
        const std::string& default_prefix, const std::string& custom_prefix);

    void operator()(amrex::MLMG& /*mlmg*/);
    void operator()(Hydro::NodalProjector& /*nodal_proj*/);
    void operator()(Hydro::MacProjector& /*mac_proj*/);

    //! Linear operator options during construction
    amrex::LPInfo& lpinfo() { return m_lpinfo; }

    // Linear operator options
    int max_order{2};

    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real rel_tol{1.0e-11};

    //! Absolute tolerance for convergence checks
    amrex::Real abs_tol{1.0e-14};

private:
    void parse_options(const std::string& /*prefix*/);

    //! Linear operator info object
    amrex::LPInfo m_lpinfo;

    // MLMG options

    //! Bottom solver type bicgstab, cg, hypre, etc.
    std::string bottom_solver_type{"bicgstab"};

    //! Options namespace for hypre solvers
    std::string hypre_namespace{"hypre"};

    //! Options of hypre interface
    std::string hypre_interface{"ij"};

    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real bottom_rel_tol{1.0e-4};

    //! Absolute tolerance for convergence checks
    amrex::Real bottom_abs_tol{-1.0};

    /** Verbosity of the MLMG and Projector interfaces
     */
    int verbose{0};
    int max_iter{200};
    int max_fmg_iters{0};
    int num_pre_smooth{2};
    int num_post_smooth{2};
    int num_final_smooth{8};
    int num_bottom_smooth{0};

    /** Set verbosity for AMReX bottom solvers
     *
     *  When using hypre as the bottom solver, this option is not passed to
     *  hypre. Instead user must set `hypre.verbose` to control the hypre solver
     *  verbosity and `hypre.bamg_verbose` to control BoomerAMG verbosity (when
     *  used as a preconditioner).
     */
    int bottom_verbose{0};

    //! Maximum iterations for the bottom solver
    int bottom_max_iter{200};

    //! If set to true, AMReX won't abort if linear solver doesn't converge
    //! within max iterations
    bool do_fixed_iters{true};

    //! option to use the nsolve bottom solver
    //! added to help with overset masking
    bool do_nsolve{false};

    //! set grid size for nsolve
#ifdef AMREX_USE_GPU
    int nsolve_grid_size{64};
#else
    int nsolve_grid_size{16};
#endif
};

} // namespace amr_wind

#endif /* MLMGOPTIONS_H */
