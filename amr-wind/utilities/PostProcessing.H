#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <memory>

#include "amr-wind/core/Factory.H"

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
    static const std::string base_identifier() { return "PostProcessBase"; }

    virtual ~PostProcessBase() = default;

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

    //! Actions to perform post regrid
    virtual void post_regrid_actions() = 0;
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
