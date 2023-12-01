#ifndef ABLUSERDEFINED_H
#define ABLUSERDEFINED_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_Gpu.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {

// signatures should match external_udf.h
using Vfun_ptr =
    void (*)(double time, double x, double y, double z, double (&ptvel)[3]);
using Tfun_ptr =
    void (*)(double time, double x, double y, double z, double& pttemp);

struct UserDefinedFunctions
{
    Vfun_ptr get_vel;
    Tfun_ptr get_temp;
};

/** Interface for ABL user-coded time-varying inlet condition
 *  \ingroup we_abl
 *
 *  This class sets velocity and temperature at mass inflow bc's
 *
 */
class ABLUserDefined
{

public:
    explicit ABLUserDefined(CFDSim& /*sim*/);

    ~ABLUserDefined();

    //! Execute initialization actions after mesh has been fully generated
    void post_init_actions();

    void pre_advance_work();

    void post_advance_work();

    void set_velocity(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        const int dcomp = 0,
        const int orig_comp = 0) const;

    void set_temperature(
        const int lev,
        const amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

private:
    const CFDSim& m_sim;
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_temperature;

    bool m_active{false};
    std::string m_inflow_lib;

#ifdef AMR_WIND_USE_EXTERNAL_UDF
    void* userfun_lib;
    UserDefinedFunctions m_udf;
#endif
};

} // namespace amr_wind

#endif /* ABLUSERDEFINED_H */
