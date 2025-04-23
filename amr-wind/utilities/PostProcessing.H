#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <memory>

#include "amr-wind/core/Factory.H"
#include "AMReX_ParmParse.H"

/**
 *  \defgroup utilities Utilities
 *  Utilities
 *
 *  This group documents various utilities available in AMR-Wind for
 *  input/output, post-processing, collecting statistics etc.
 */

namespace amr_wind {

class CFDSim;

/** Abstract representation of a post-processing utility
 *  \ingroup utilities
 *
 *  Provides an interface for a runtime-selectable post-processing action (e.g.,
 *  data sampling, volume/surface integration, etc.) with the main solver.
 *
 *  All post-processing utilities must derive from this class.
 */
class PostProcessBase
    : public Factory<PostProcessBase, CFDSim&, const std::string&>
{
public:
    static std::string base_identifier() { return "PostProcessBase"; }

    ~PostProcessBase() override = default;

    virtual void pre_init_actions() = 0;

    /** Initialize the post-processing utility
     *
     *  When this method is called, the mesh and the solution fields are fully
     *  initialized
     */
    virtual void initialize() = 0;

    /** Perform actions at the end of each timestep
     */
    virtual void post_advance_work() = 0;

    /** Perform output routine
     */
    virtual void output_actions() = 0;

    //! Actions to perform post regrid
    virtual void post_regrid_actions() = 0;

    void populate_output_parameters(amrex::ParmParse& pp)
    {
        pp.query("output_interval", m_out_interval);
        pp.query("output_time_interval", m_out_time_interval);
        pp.query("output_delay", m_out_delay);
        pp.query("output_time_delay", m_out_time_delay);
        pp.query("enforce_output_time_dt", m_enforce_dt);
        pp.query("enforce_output_time_dt_reltol", m_enforce_dt_tol);

        if (!pp.contains("output_interval")) {
            pp.query("output_frequency", m_out_interval);
            if (m_out_interval > 0) {
                amrex::Print()
                    << "WARNING: legacy input argument of output_frequency has "
                       "been specified; the code will function as normal, but "
                       "the new input argument of output_interval is "
                       "preferred.\n";
            }
        }

        if (m_out_interval > 0 && m_out_time_interval > 0.0) {
            amrex::Abort(
                "output_interval (or output_frequency) and "
                "output_time_interval are both specified. timestep- and "
                "time-based output parameters should not be used together; "
                "please only specify one.");
        }

        if (m_out_interval <= 0 && m_out_time_interval <= 0.0) {
            amrex::Abort(
                "output_interval (or output_frequency) and "
                "output_time_interval are both not specified. please specify "
                "one of these input arguments so that this postprocessing tool "
                "can output data.");
        }

        if (m_out_time_interval <= 0.0 && m_enforce_dt) {
            amrex::Abort(
                "enforce_output_time_dt is true, but output_time_interval "
                "has not been specified.");
        }
    }

    //! Determine if output will happen this time step
    bool do_output_now(
        const int ntime,
        const amrex::Real time,
        const amrex::Real dt,
        amrex::Real tol) const
    {

        // Choose smallest tolerance between incoming and local
        tol = std::min(tol, m_out_time_tol * dt);

        const bool step_based =
            (m_out_interval > 0) &&
            ((ntime - m_out_delay >= 0) && (ntime % m_out_interval == 0));
        const bool time_based =
            (m_out_time_interval > 0.) &&
            ((time + tol - m_out_time_delay >= 0.) &&
             (time + tol -
                  std::floor((time + tol) / m_out_time_interval) *
                      m_out_time_interval <
              dt));

        return (step_based || time_based);
    }

    amrex::Real output_time_interval() const { return m_out_time_interval; }
    amrex::Real output_time_delay() const { return m_out_time_delay; }
    amrex::Real enforce_dt_tolerance() const { return m_enforce_dt_tol; }
    bool enforce_dt() const { return m_enforce_dt; }

protected:
    //! Time step interval for output
    int m_out_interval{-1};
    //! Time interval for output
    amrex::Real m_out_time_interval{-1.};
    //! Time step delay for output
    int m_out_delay{0};
    //! Time delay for output
    amrex::Real m_out_time_delay{0.};
    //! Relative (to dt) tolerance for time interval output
    amrex::Real m_out_time_tol{1e-8};
    //! Relative (to m_out_time_interval) tolerance for enforcing dt
    amrex::Real m_enforce_dt_tol{1e-3};
    //! Flag for enforcing time step based on output time
    bool m_enforce_dt{false};
};

/** A collection of post-processing instances that are active during a
 *  simulation
 *
 *  Upon initialization, the class reads the ``incflo.post_processing``
 *  parameter to determine the utilities that must be activated and uses runtime
 *  selection to initialize them
 *
 *  \ingroup utilities
 */
class PostProcessManager
{
public:
    explicit PostProcessManager(CFDSim& sim);

    ~PostProcessManager() = default;

    void pre_init_actions();

    /** Initialize post-processing utilities
     *
     *  Reads user inputs and instantiates all the registered post-processing
     *  utilities. Executes the registered utilities once before starting the
     *  time integration
     */
    void post_init_actions();

    //! Call all registered utilities to perform actions after a timestep
    void post_advance_work();

    void post_regrid_actions();

private:
    CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<PostProcessBase>> m_post;
};

} // namespace amr_wind

#endif /* POSTPROCESSING_H */
