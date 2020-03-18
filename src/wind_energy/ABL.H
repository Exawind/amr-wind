#ifndef ABL_H
#define ABL_H

#include <memory>

#include "Physics.H"
#include "ABLFieldInit.H"

namespace amr_wind {

/** Atmospheric Boundary Layer physics
 */
class ABL : public Physics
{
public:
    ABL();

    ~ABL() = default;

    std::string physics_name() const override { return "ABL"; }

    void initialize_fields(
        const amrex::Geometry& geom,
        amrex::MultiFab& density,
        amrex::MultiFab& velocity,
        amrex::MultiFab& pressure,
        amrex::MultiFab& scalars) const override;

    const ABLFieldInit& abl_initializer() const { return *m_field_init; }

private:
    //! ABL field initializer instance
    std::unique_ptr<ABLFieldInit> m_field_init;

    //! True if Boussinesq buoyancy term is active
    bool m_has_boussinesq{true};

    //! True if Coriolis forcing term is active
    bool m_has_coriolis{true};

    //! True if driving pressure gradient term is activate
    bool m_has_driving_dpdx{true};
};

}

#endif /* ABL_H */
