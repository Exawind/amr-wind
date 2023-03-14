#ifndef WING_OPS_H
#define WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"

namespace amr_wind::actuator {
namespace wing {

/** Parse common data for modeling fixed wings as actuators
 */
void read_inputs(
    WingBaseData& /*wdata*/, ActInfo& /*info*/, const utils::ActParser& /*pp*/);

/** Initialize core data structures when modeling fixed wings
 */
void init_data_structures(WingBaseData& /*wdata*/, ActGrid& /*grid*/);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/,
    const amrex::Real /*time*/);

template <typename T>
ComponentView make_component_view(typename T::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();
    ComponentView view;
    view.pos = ::amr_wind::utils::slice(grid.pos, 0, meta.num_pts);
    view.vel_pos = ::amr_wind::utils::slice(grid.vel_pos, 0, meta.num_pts);
    view.force = ::amr_wind::utils::slice(grid.force, 0, meta.num_pts);
    view.epsilon = ::amr_wind::utils::slice(grid.epsilon, 0, meta.num_pts);
    view.orientation =
        ::amr_wind::utils::slice(grid.orientation, 0, meta.num_pts);
    view.chord = ::amr_wind::utils::slice(meta.chord, 0, meta.num_pts);
    view.vel_rel = ::amr_wind::utils::slice(meta.vel_rel, 0, meta.num_pts);
    view.vel = ::amr_wind::utils::slice(grid.vel, 0, meta.num_pts);

    return view;
}

} // namespace wing

namespace ops {

template <typename T>
const typename AirfoilTraits<T>::AirfoilLookup&
airfoil_lookup(typename T::DataType& data);

template <typename ActTrait, typename SrcTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_pos>::type>
{
    void operator()(typename ActTrait::DataType& /*unused*/) {}
};

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_vel>::type>
{

    void operator()(typename ActTrait::DataType& data)
    {
        auto& meta = data.meta();
        if (meta.fllc) {
            FLLCOp()(meta.component_view, *(meta.fllc));
        }
    }
};

template <typename ActTrait>
struct ComputeForceOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_force>::type>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        const auto& chord = wdata.chord;
        const auto& aflookup = airfoil_lookup<ActTrait>(data);

        amrex::Real total_lift = 0.0;
        amrex::Real total_drag = 0.0;
        for (int ip = 0; ip < npts; ++ip) {
//            const auto& tmat = grid.orientation[ip];
            // Effective velocity at the wing control point in local frame
//            auto wvel = tmat & grid.vel[ip];
            // Set spanwise component to zero to get a pure 2D velocity
//            wvel.y() = 0.0;

            // Build the local reference frame
            vs::Vector wspan = wdata.end - wdata.start;
            wspan = wspan.unit();
  
            // Use the global coord to orient the blade
            // This works only for inflow in the x direction
            auto blade_x = wdata.blade_x;
            auto blade_y = wspan;
            auto blade_z = (blade_x^blade_y).unit();

//std::cout << "blade_x" << std::endl;
//std::cout << blade_x[0] << std::endl;
//std::cout << blade_x[1] << std::endl;
//std::cout << blade_x[2] << std::endl;
//std::cout << "blade_y" << std::endl;
//std::cout << blade_y[0] << std::endl;
//std::cout << blade_y[1] << std::endl;
//std::cout << blade_y[2] << std::endl;
//std::cout << "blade_z" << std::endl;
//std::cout << blade_z[0] << std::endl;
//std::cout << blade_z[1] << std::endl;
//std::cout << blade_z[2] << std::endl;

// Quaternion used to rotate (x, y, z) onto blade local reference frame
//auto uv = amr_wind::utils::degrees(1. + (local_y & wspan));
//auto tmat2 = vs::quaternion(blade_y,  uv);
// Rotate the other axes 
//auto blade_x = tmat2 & local_x;
//auto blade_z = tmat2 & local_z;



//            const auto vmag = vs::mag(wvel);
//            const auto aoa = std::atan2(wvel.z(), wvel.x());

vs::Vector windvector;
windvector[0] = grid.vel[ip] & blade_x;
windvector[1] = 0; //grid.vel[ip] & blade_y;
windvector[2] = grid.vel[ip] & blade_z;

const auto vmag = vs::mag(windvector);
//const auto aoa = amr_wind::utils::degrees(std::atan2(windvector[2], windvector[0])) + wdata.pitch;
const auto aoa = std::atan2(windvector[2], windvector[0]) + amr_wind::utils::radians(wdata.pitch);
//std::cout << aoa << std::endl;
//std::cout << "Pitch " << wdata.pitch << std::endl;

            // Make up some Cl, Cd values
            amrex::Real cl, cd;
            aflookup(aoa, cl, cd); // This function calls aoa in radians (-pi to pi) 

            // Assume unit chord
            const auto qval = 0.5 * vmag * vmag * chord[ip] * dx[ip];
            const auto lift = qval * cl;
            const auto drag = qval * cd;
            // Determine unit vector parallel and perpendicular to velocity
            // vector
//            const auto drag_dir = wvel.unit() & tmat;
//            const auto lift_dir = drag_dir ^ tmat.y();

// Directions 
const auto drag_dir = (blade_x * windvector.x() + blade_z * windvector.z() ).unit();
const auto lift_dir = (drag_dir ^ blade_y).unit();

            // Compute force on fluid from this section of wing
            grid.force[ip] = -(lift_dir * lift + drag * drag_dir);

            // Assign values for output
            wdata.vel_rel[ip] = windvector;
            wdata.aoa[ip] = amr_wind::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;

        if (wdata.fllc) {
            if (!(wdata.fllc->initialized)) {
                wdata.component_view =
                    amr_wind::actuator::wing::make_component_view<ActTrait>(
                        data);
                FLLCInit(
                    *(wdata.fllc.get()), wdata.component_view,
                    wdata.epsilon_chord[0]);
            }
        }
    }
};

template <typename ActTrait, typename SrcTrait>
struct ProcessOutputsOp<
    ActTrait,
    SrcTrait,
    typename std::enable_if<
        std::is_base_of<WingType, ActTrait>::value &&
        UseDefaultOp<ActTrait, ActSrcLine>::process_outputs>::type>
{
private:
    typename ActTrait::DataType& m_data;

    std::string m_out_dir;
    std::string m_nc_filename;
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename ActTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        wing::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        wing::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace amr_wind::actuator

#endif /* WING_OPS_H */
