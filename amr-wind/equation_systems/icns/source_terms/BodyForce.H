#ifndef BODYFORCE_H
#define BODYFORCE_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind {
namespace pde {
namespace icns {

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
    const Field& m_density;
   
    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_body_force{{0.0, 0.0, 0.0}};

    //! Body Force Type
    std::string m_type;

    //! Angular frequency used in the oscillatory forcing
    amrex::Real m_omega;

    //! Time
    const SimTime& m_time;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* BODYFORCE_H */
