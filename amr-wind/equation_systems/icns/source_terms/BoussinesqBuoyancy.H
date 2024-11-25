#ifndef BOUSSINESQBUOYANCY_H
#define BOUSSINESQBUOYANCY_H

#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term
 *  \ingroup icns_src we_abl
 *
 *  \f[
 *    S = \beta g \left( T_\mathrm{ref} - T \right)
 *  \f]
 */
class BoussinesqBuoyancy : public MomentumSource::Register<BoussinesqBuoyancy>
{
public:
    static std::string identifier() { return "BoussinesqBuoyancy"; }

    explicit BoussinesqBuoyancy(const CFDSim& sim);

    ~BoussinesqBuoyancy() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_temperature;
    const Field* m_vof;

    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};

    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};

    //! Check for VOF
    bool m_is_vof{false};

  //! Transport model
  std::unique_ptr<transport::TransportModel> m_transport;
};
} // namespace amr_wind::pde::icns

#endif /* BOUSSINESQBUOYANCY_H */
