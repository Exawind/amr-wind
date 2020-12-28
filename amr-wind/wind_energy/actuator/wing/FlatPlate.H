#ifndef FLATPLATE_H
#define FLATPLATE_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {
namespace actuator {
namespace wing {

struct ThinAirfoil : public AirfoilLookup
{
    amrex::Real cd_fac{0.0};

    ThinAirfoil() = default;
    ThinAirfoil(amrex::Real cd_fac_in) : cd_fac(cd_fac_in) {}

    AMREX_FORCE_INLINE
    void operator()(
        const amrex::Real aoa, amrex::Real& cl, amrex::Real& cd) const override
    {
        cl = ::amr_wind::utils::two_pi() * aoa;
        cd = cd_fac * cl * std::sin(aoa);
    }
};

} // namespace wing

struct FlatPlateData : public WingBaseData
{
    amrex::Real chord_inp{1.0};
    amrex::Real cd_factor{0.0};

    wing::ThinAirfoil aflookup;
};

struct FlatPlate : public WingType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = FlatPlateData;
    using DataType = ActDataHolder<FlatPlate>;

    static const std::string identifier() { return "FlatPlate"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* FLATPLATE_H */
