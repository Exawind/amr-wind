#ifndef ABLANELASTIC_H
#define ABLANELASTIC_H
#include "amr-wind/CFDSim.H"

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
    void initialize_data();
    bool is_anelastic() const { return m_is_anelastic; }

private:
    const amrex::AmrCore& m_mesh;

    bool m_is_anelastic{false};
};
} // namespace amr_wind
#endif
