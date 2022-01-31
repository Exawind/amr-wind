#ifndef COMPLEX_TERRAIN_OPS_H
#define COMPLEX_TERRAIN_OPS_H

#include "amr-wind/immersed_boundary/complex_terrain/ComplexTerrain.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/fvm/gradient.H"

#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace ib {
namespace complex_terrain {

/** Parse common data for modeling non-blade bluff bodies via an immersed
 * bounary
 */
void read_inputs(
    ComplexTerrainBaseData&, IBInfo&, const ::amr_wind::utils::MultiParser&);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(CFDSim&);

/** Set the velocity inside the IB Complex Terrain
 */
void apply_dirichlet_vel(CFDSim&, const amrex::Vector<amrex::Real>& vel_bc);

void prepare_netcdf_file(
    const std::string&, const ComplexTerrainBaseData&, const IBInfo&);

void write_netcdf(
    const std::string&,
    const ComplexTerrainBaseData&,
    const IBInfo&,
    const amrex::Real);

} // namespace complex_terrain

namespace ops {

template <typename GeomTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename GeomTrait>
struct UpdatePosOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<ComplexTerrainType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct UpdateVelOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<ComplexTerrainType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ib::IB::update_velocities");

        const auto& wdata = data.meta();
        auto& sim = data.sim();

        complex_terrain::apply_dirichlet_vel(sim, wdata.vel_bc);
    }
};

template <typename GeomTrait>
struct ComputeForceOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<ComplexTerrainType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct ProcessOutputsOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<ComplexTerrainType, GeomTrait>::value>::type>
{
private:
    typename GeomTrait::DataType& m_data;

    std::string m_out_dir{""};
    std::string m_nc_filename{""};
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename GeomTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const ::amr_wind::utils::MultiParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        complex_terrain::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        complex_terrain::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), time.new_time());
    }
};

} // namespace ops

} // namespace ib
} // namespace amr_wind

#endif /* COMPLEX_TERRAIN_OPS_H */
