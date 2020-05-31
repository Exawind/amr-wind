#ifndef CORIOLISFORCING_H
#define CORIOLISFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind {
namespace pde {
namespace icns {

class CoriolisForcing : public MomentumSource::Register<CoriolisForcing>
{
public:
    static const std::string identifier() { return "CoriolisForcing"; }

    CoriolisForcing(const CFDSim&);

    virtual ~CoriolisForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_velocity;

    ///@{
    /** Orientation of the grid w.r.t. to planetary coordinate system
     */
    amrex::Vector<amrex::Real> m_east{{1.0, 0.0, 0.0}};
    amrex::Vector<amrex::Real> m_north{{0.0, 1.0, 0.0}};
    amrex::Vector<amrex::Real> m_up{{1.0, 0.0, 1.0}};
    ///@}

    //! Latitude where the Coriolis forcing is computed (+ve northern
    //! hemisphere)
    amrex::Real m_latitude{0.0};

    amrex::Real m_sinphi{0.0};
    amrex::Real m_cosphi{0.0};

    //! `2.0 * \Omega`
    amrex::Real m_coriolis_factor{0.0};
};

} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* CORIOLISFORCING_H */
