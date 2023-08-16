#ifndef FIXEDWING_H
#define FIXEDWING_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"

namespace amr_wind::actuator {

struct FixedWingData : public WingBaseData
{
    RealList span_locs{0.0, 1.0};
    RealList chord_inp{1.0, 1.0};
    std::string airfoil_file;
    std::string airfoil_type{"openfast"};

    std::unique_ptr<AirfoilTable> aflookup;
};

struct FixedWing : public WingType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = FixedWingData;
    using DataType = ActDataHolder<FixedWing>;

    static std::string identifier() { return "FixedWing"; }
};

} // namespace amr_wind::actuator

#endif /* FIXEDWING_H */
