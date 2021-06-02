#ifndef TURBINE_TYPES_H
#define TURBINE_TYPES_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"

namespace amr_wind {
namespace actuator {

struct TurbineInfo : public ActInfo
{
    //! Base location of this turbine
    vs::Vector base_pos;

    //! Rotor diameter for this turbine
    amrex::Real rotor_diameter{0.0};

    //! Hub height of this turbine
    amrex::Real hub_height{0.0};

    TurbineInfo(const std::string& label_in, const int id_in)
        : ActInfo(label_in, id_in)
    {}
};

struct ComponentView
{
    VecSlice pos;
    VecSlice force;
    VecSlice epsilon;
    TensorSlice orientation;

    VecSlice vel_pos;
    VecSlice vel;

    RealSlice chord;
};

using DeviceVecComponent = amrex::Gpu::DeviceVector<ComponentView>;
struct TurbineBaseData
{
    //! Number of blades
    int num_blades{3};

    //! Number of actuator nodes per blade
    int num_pts_blade;

    //! Number of actuator nodes for tower
    int num_pts_tower;

    RealList chord;

    //! Center of rotation for the rotor
    vs::Vector rot_center;

    /** Reference frame for the rotor plane
     *
     *  x - Pointing downwind along shaft/rotation axis
     *  y - lateral direction
     *  z - normal direction (pointing mostly along z-direction)
     */
    vs::Tensor rotor_frame;

    //! General epsilon provided
    vs::Vector eps_inp;

    //! Gaussian smearing factor for the blade
    vs::Vector eps_chord;

    //! Minimum epsilon when using chord based epsilon
    vs::Vector eps_min;

    //! Gaussian smearing factor for tower
    vs::Vector eps_tower;

    //! Drag coefficient for nacelle
    amrex::Real nacelle_cd{0.0};

    //! Wetted surface area for nacelle
    amrex::Real nacelle_area{0.0};

    std::vector<ComponentView> blades;
    ComponentView tower;
    ComponentView hub;
};

struct TurbineType : public ActuatorType
{};

} // namespace actuator
} // namespace amr_wind

#endif /* TURBINE_TYPES_H */
