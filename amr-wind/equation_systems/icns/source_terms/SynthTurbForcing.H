#ifndef SYNTHTURBFORCING_H
#define SYNTHTURBFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Synthetic turbulence momentum source term
 *  \ingroup icns_src synth_turb
 */
class SynthTurbForcing : public MomentumSource::Register<SynthTurbForcing>
{
public:
    static const std::string identifier() { return "SynthTurbForcing"; }

    explicit SynthTurbForcing(const CFDSim&);

    virtual ~SynthTurbForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_turb_force;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* SYNTHTURBFORCING_H */
