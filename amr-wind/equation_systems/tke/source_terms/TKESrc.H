#ifndef TKESRC_H
#define TKESRC_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind {
namespace pde {
namespace tke {

/** Source terms for TKE equation
 *   \ingroup turb_model
 *
 *  Works for all TKE equations with source and dissipation terms
 */
class TKESrc : public TKESource::Register<TKESrc>
{
public:
    static const std::string identifier() { return "TKESrc"; }

    TKESrc(const CFDSim&);

    virtual ~TKESrc();

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

#endif /* TKESRC_H */
