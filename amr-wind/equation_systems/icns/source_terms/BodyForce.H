#ifndef BODYFORCE_H
#define BODYFORCE_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Apply a user-defined body force as a source term
 *  \ingroup icns_src
 *
 *  The implementation can apply either a constant source term, or a sinusoidal
 *  term as a function of time.
 */
class BodyForce : public MomentumSource::Register<BodyForce>
{
public:
    static const std::string identifier() { return "BodyForce"; }

    BodyForce(const CFDSim&);

    virtual ~BodyForce();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    //! Time
    const SimTime& m_time;

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_body_force{{0.0, 0.0, 0.0}};

    //! Body Force Type
    std::string m_type{"constant"};

    //! Angular frequency used in the oscillatory forcing
    amrex::Real m_omega{0.0};
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* BODYFORCE_H */
