#ifndef WALLFUNCTION_H
#define WALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {

/** Wall-function computations for LES simulations
 *  \ingroup utilities
 *
 *  This class performs the necessary computations at the beginning of
 *  predictor/corrector steps. The actual BC population in ghost cells is
 *  performed by VelWallFunc BC interface class.
 */
class WallFunction
{
public:
    explicit WallFunction(const CFDSim& sim);

    amrex::Real utau() const { return m_utau; }

    ~WallFunction() = default;

private:
    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    amrex::Real m_utau{0.0};
};

/** Applies a shear-stress value at the domain boundary
 *  \ingroup field_bc utilities
 *
 *  \sa WallFunction
 */
class VelWallFunc : public FieldBCIface
{
public:
    VelWallFunc(Field& velocity, const WallFunction& wall_func);

    void operator()(Field& velocity, const FieldState rho_state) override;

    static void wall_model(
        Field& velocity, const FieldState rho_state, const amrex::Real utau);

private:
    const WallFunction& m_wall_func;
    std::string m_wall_shear_stress_type{"constant"};
};
} // namespace amr_wind

#endif /* WALLFUNCTION_H */
