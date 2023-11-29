#ifndef ABLCODEDINLET_H
#define ABLCODEDINLET_H

#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_Gpu.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind {

// signatures should match coded_inlet_fields.h
using Vfun_ptr =
    void (*)(double time, double x, double y, double z, double (&ptvel)[3]);
using Tfun_ptr =
    void (*)(double time, double x, double y, double z, double& pttemp);

struct UserCoded
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
class ABLCodedInlet
{

public:
    explicit ABLCodedInlet(CFDSim& /*sim*/);

    ~ABLCodedInlet();

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
    std::string m_user_lib;

    void* userfun_lib;
    UserCoded m_user_coded;
};

} // namespace amr_wind

#endif /* ABLCODEDINLET_H */
