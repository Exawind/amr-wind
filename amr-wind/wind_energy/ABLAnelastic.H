#ifndef ABLANELASTIC_H
#define ABLANELASTIC_H
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

/** Interface for ABL anelastic
 *  \ingroup we_abl
 *
 *  This class performs the necessary operations for anelastic ABL
 *
 */
class ABLAnelastic
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    explicit ABLAnelastic(CFDSim& /*sim*/);

    void post_init_actions();
    void post_regrid_actions();
    void initialize_data();
    bool is_anelastic() const { return m_is_anelastic; }

private:
    const CFDSim& m_sim;

    bool m_is_anelastic{false};

    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    amrex::Real m_rho0_const{1.0};

    amrex::Real m_atmospheric_pressure{1.01325e5};

    const int m_axis{2};
};
} // namespace amr_wind
#endif
