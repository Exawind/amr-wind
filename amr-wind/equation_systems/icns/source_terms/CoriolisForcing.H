#ifndef CORIOLISFORCING_H
#define CORIOLISFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind::pde::icns {

/** Account for Coriolis forcing due to earth's rotation
 *
 *
 * Coriolis forcing source term is the negative of the cross
 * product of the Coriolis factor with the velocity field vector
 * Rotating the coriolis factor into the local frame of reference
 * from the earths frame of reference results in two components
 * coriolis factor(cosphi, sinphi) in the y and z directions respectively.
 * Taking the cross product of coriolis factor(cosphi, sinphi) with
 * the velocity field vector u(x,y,z) results in three source term components
 * src(x,y,z).
 * The vertical component of the Coriolis source term is typically neglected as
 * it is unphysical and is also, negligble compared to the other vertical
 * component source terms.
 *
 * There is an option to include the vertical Coriolis
 * source term by turning on three_ComponentForcing = true.
 * three_ComponentForcing is by default set to false.
 * This removes the contributions of the Coriolis factor(cosphi) component
 * from the Coriolis source term components.
 *
 *
 *  \ingroup icns_src we_abl
 *
 */
class CoriolisForcing : public MomentumSource::Register<CoriolisForcing>
{
public:
    static std::string identifier() { return "CoriolisForcing"; }

    explicit CoriolisForcing(const CFDSim& /*sim*/);

    ~CoriolisForcing() override;

    void operator()(
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
    amrex::Vector<amrex::Real> m_up{{0.0, 0.0, 1.0}};
    ///@}

    //! Latitude where the Coriolis forcing is computed (+ve northern
    //! hemisphere)
    amrex::Real m_latitude{0.0};

    amrex::Real m_sinphi{0.0};
    amrex::Real m_cosphi{0.0};

    //! Default three component forcing
    bool m_three_dimensional_forcing{false};

    //! `2.0 * \Omega`
    amrex::Real m_coriolis_factor{0.0};
};

} // namespace amr_wind::pde::icns

#endif /* CORIOLISFORCING_H */
