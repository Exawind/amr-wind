#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::temperature {

/** Adds the non-linear terms from the Kosovic Model to Temperature Equation..
 *
 *  \ingroup temperature_src NonLinearSGS
 *
 *  \sa NonLinearSGS
 */
class DragTempForcing : public TemperatureSource::Register<DragTempForcing>
{
public:
    static std::string identifier() { return "DragTempForcing"; }

    explicit DragTempForcing(const CFDSim& sim);

    ~DragTempForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const amrex::AmrCore& m_mesh;
    Field& m_temperature;
    Field& m_terrainBlank;
    amrex::Real m_drag{10.0};
};

} // namespace amr_wind::pde::temperature

#endif /* NONLINEARSGSTERM_H */
