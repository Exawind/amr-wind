#ifndef Krans_H
#define Krans_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind::pde::tke {

/** TKE source term based on Axell 2011 paper
 *  \ingroup tke_src turb_model we_abl
 */
class Krans : public TKESource::Register<Krans>
{
public:
    static std::string identifier() { return "Krans"; }

    explicit Krans(const CFDSim& /*sim*/);

    ~Krans() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_tke;
  amrex::Real m_heat_flux{0.0};
  amrex::Real m_ref_temp{300.0};
  const SimTime& m_time;
  const CFDSim& m_sim;
  const amrex::AmrCore& m_mesh;
  const Field& m_velocity;
};


} // namespace amr_wind::pde::tke

#endif /* Krans_H */
