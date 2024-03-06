#ifndef NONLINEARSGSTERM_H
#define NONLINEARSGSTERM_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde::icns {

/** Body forces introduced by turbines modeled as actuators in flow field.
 *
 *  \ingroup icns_src actuator
 *
 *  \sa Acuator
 */
class NonLinearSGSTerm : public MomentumSource::Register<NonLinearSGSTerm>
{
public:
    static std::string identifier() { return "NonLinearSGSTerm"; }

    explicit NonLinearSGSTerm(const CFDSim& sim);

    ~NonLinearSGSTerm() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_divNij;
};

} // namespace amr_wind::pde::icns

#endif /* ACTUATORFORCING_H */
