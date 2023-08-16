#include "amr-wind/core/MLMGOptions.H"

#include "hydro_MacProjector.H"
#include "hydro_NodalProjector.H"

namespace amr_wind {

MLMGOptions::MLMGOptions(const std::string& prefix) { parse_options(prefix); }

MLMGOptions::MLMGOptions(
    const std::string& default_prefix, const std::string& custom_prefix)
{
    parse_options(default_prefix);
    parse_options(custom_prefix);
}

void MLMGOptions::parse_options(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);

    // LPInfo options
    pp.query("do_agglomeration", m_lpinfo.do_agglomeration);
    pp.query("do_consolidation", m_lpinfo.do_consolidation);
    pp.query("do_semicoarsening", m_lpinfo.do_semicoarsening);
    pp.query("agg_grid_size", m_lpinfo.agg_grid_size);
    pp.query("con_grid_size", m_lpinfo.con_grid_size);
    pp.query("max_coarsening_level", m_lpinfo.max_coarsening_level);
    pp.query("max_semicoarsening_level", m_lpinfo.max_semicoarsening_level);

    // LinOp options
    pp.query("max_order", max_order);

    // MLMG options
    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);

    pp.query("maxiter", max_iter);
    pp.query("mg_rtol", rel_tol);
    pp.query("mg_atol", abs_tol);
    pp.query("fmg_maxiter", max_fmg_iters);
    pp.query("num_pre_smooth", num_pre_smooth);
    pp.query("num_post_smooth", num_post_smooth);
    pp.query("num_final_smooth", num_final_smooth);
    pp.query("num_bottom_smooth", num_bottom_smooth);

    pp.query("do_fixed_iters", do_fixed_iters);

    pp.query("bottom_maxiter", bottom_max_iter);
    pp.query("bottom_rtol", bottom_rel_tol);
    pp.query("bottom_atol", bottom_abs_tol);
    pp.query("bottom_solver", bottom_solver_type);
    pp.query("hypre_namespace", hypre_namespace);
    pp.query("hypre_interface", hypre_interface);
    pp.query("do_nsolve", do_nsolve);
    pp.query("nsolve_grid_size", nsolve_grid_size);
}

void MLMGOptions::operator()(amrex::MLMG& mlmg)
{
    mlmg.setVerbose(verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iters);

    if (do_fixed_iters) {
        mlmg.setFixedIter(max_iter);
    }

    mlmg.setNSolve(static_cast<int>(do_nsolve));
    mlmg.setNSolveGridSize(nsolve_grid_size);
    mlmg.setPreSmooth(num_pre_smooth);
    mlmg.setPostSmooth(num_post_smooth);
    mlmg.setFinalSmooth(num_final_smooth);
    mlmg.setBottomSmooth(num_bottom_smooth);

    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setBottomTolerance(bottom_rel_tol);
    mlmg.setBottomToleranceAbs(bottom_abs_tol);

    if (bottom_solver_type == "smoother") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
    } else if (bottom_solver_type == "bicg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::bicgstab);
    } else if (bottom_solver_type == "cg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::cg);
    } else if (bottom_solver_type == "bicgcg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::bicgcg);
    } else if (bottom_solver_type == "cgbicg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::cgbicg);
    } else if (bottom_solver_type == "hypre") {
#ifdef AMREX_USE_HYPRE
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::hypre);

        mlmg.setHypreOptionsNamespace(hypre_namespace);
        if (hypre_interface == "ij")
            mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
        else if (hypre_interface == "semi_structured")
            mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
        else if (hypre_interface == "structured")
            mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
        else
            amrex::Abort(
                "Invalid hypre interface. Valid options: ij semi_structured "
                "structured");
#else
        amrex::Abort("AMR-Wind was not built with hypre support");
#endif
    }
}

void MLMGOptions::operator()(Hydro::MacProjector& mac_proj)
{
    operator()(mac_proj.getMLMG());
}

void MLMGOptions::operator()(Hydro::NodalProjector& nodal_proj)
{
    // Only IJ interface supported for NodalProjector
    hypre_interface = "ij";

    if ((m_lpinfo.max_coarsening_level > -1) &&
        (bottom_solver_type == "hypre")) {
        nodal_proj.getLinOp().setCoarseningStrategy(
            amrex::MLNodeLinOp::CoarseningStrategy::RAP);
    }

    nodal_proj.setVerbose(verbose);
    operator()(nodal_proj.getMLMG());
}

} // namespace amr_wind
