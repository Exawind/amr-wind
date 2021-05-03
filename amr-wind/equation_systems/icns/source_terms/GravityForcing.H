#ifndef GravityForcing_H
#define GravityForcing_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Gravity Forcing term
 *  \ingroup icns_src
 */
class GravityForcing : public MomentumSource::Register<GravityForcing>
{
public:
    static const std::string identifier() { return "GravityForcing"; }

    explicit GravityForcing(const CFDSim& sim);

    virtual ~GravityForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* GravityForcing_H */
