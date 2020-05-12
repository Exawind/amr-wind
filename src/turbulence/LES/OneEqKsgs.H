#ifndef ONEEQKSGS_H
#define ONEEQKSGS_H

#include <string>
#include "TurbModelBase.H"

namespace amr_wind {
namespace turbulence {

/** Base class for 1-Equation subgrid scale TKE turbulence model
 */
template <typename Transport>
class OneEqKsgs : public TurbModelBase<Transport>
{
public:
    OneEqKsgs(CFDSim& sim);

    virtual ~OneEqKsgs();

protected:
    Field& m_vel;

    //! Turbulence constant
    amrex::Real m_Ceps{0.93};
};

template <typename Transport>
class OneEqKsgsM84 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsM84-" + Transport::identifier();
    }

    OneEqKsgsM84(CFDSim& sim);

    virtual ~OneEqKsgsM84();

    virtual std::string model_name() const override { return "OneEqKsgsM84"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;
};

template <typename Transport>
class OneEqKsgsS94 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsS94-" + Transport::identifier();
    }

    OneEqKsgsS94(CFDSim& sim);

    virtual ~OneEqKsgsS94();

    virtual std::string model_name() const override { return "OneEqKsgsS94"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;
};

} // namespace turbulence
} // namespace amr_wind

#endif /* ONEEQKSGS_H */
