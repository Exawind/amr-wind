#ifndef KWSSTSRC_H
#define KWSSTSRC_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind {
namespace pde {
namespace tke {

/** Source terms for TKE equation in K-w-SST turbulence model
 *   \ingroup turb_model
 *
 *  Works for all TKE equations with production and dissipation terms
 */
class KwSSTSrc : public TKESource::Register<KwSSTSrc>
{
public:
    static const std::string identifier() { return "KwSSTSrc"; }

    explicit KwSSTSrc(const CFDSim&);

    virtual ~KwSSTSrc();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_shear_prod;
    Field& m_diss;
};

} // namespace tke
} // namespace pde
} // namespace amr_wind

#endif /* KWSSTSRC_H */
