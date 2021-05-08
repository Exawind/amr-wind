#ifndef BLUFF_BODY_OPS_H
#define BLUFF_BODY_OPS_H

#include "amr-wind/immersed_boundary/geometry/BluffBody.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IBParser.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

/** Parse common data for modeling non-blade bluff bodies via an immersed
 * bounary
 */
void read_inputs(BluffBodyBaseData&, IBInfo&, const utils::IBParser&);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(BluffBodyBaseData&, IBGrid&);

void prepare_netcdf_file(
    const std::string&, const BluffBodyBaseData&, const IBInfo&, const IBGrid&);

void write_netcdf(
    const std::string&,
    const BluffBodyBaseData&,
    const IBInfo&,
    const IBGrid&,
    const amrex::Real);

} // namespace bluff_body

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
        std::is_base_of<BluffBodyType, GeomTrait>::value &&
        UseDefaultOp<GeomTrait>::update_pos>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct UpdateVelOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value &&
        UseDefaultOp<GeomTrait>::update_vel>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct ComputeForceOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value &&
        UseDefaultOp<GeomTrait>::compute_force>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct ProcessOutputsOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value &&
        UseDefaultOp<GeomTrait>::process_outputs>::type>
{
private:
    typename GeomTrait::DataType& m_data;

    std::string m_out_dir{""};
    std::string m_nc_filename{""};
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename GeomTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const utils::IBParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        bluff_body::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        bluff_body::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* BLUFF_BODY_OPS_H */
