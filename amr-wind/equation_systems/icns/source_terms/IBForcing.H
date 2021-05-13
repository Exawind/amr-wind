#ifndef IBFORCING_H
#define IBFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Body forces introduced by immersed boundaries a.
 *
 *  \ingroup icns_src immersed_boundary
 *
 *  \sa IB
 */
class IBForcing : public MomentumSource::Register<IBForcing>
{
public:
    static std::string identifier() { return "IBForcing"; }

    explicit IBForcing(const CFDSim& sim);

    virtual ~IBForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_ib_src;
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* ACTUATORFORCING_H */
