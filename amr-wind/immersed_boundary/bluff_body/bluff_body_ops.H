#ifndef BLUFF_BODY_OPS_H
#define BLUFF_BODY_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/core/MultiParser.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

/** Parse common data for modeling non-blade bluff bodies via an immersed
 * bounary
 */
void read_inputs(
    BluffBodyBaseData&, IBInfo&, const ::amr_wind::utils::MultiParser&);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(BluffBodyBaseData&);

void prepare_netcdf_file(
    const std::string&, const BluffBodyBaseData&, const IBInfo&);

void write_netcdf(
    const std::string&,
    const BluffBodyBaseData&,
    const IBInfo&,
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
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct UpdateVelOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ib::IB::update_velocities");

        auto& wdata = data.meta();
        auto& sim = data.sim();

        if (wdata.has_wall_model) {
            amrex::Abort("IB wall modelling is not available yet.");
        }

        const int nlevels = sim.repo().num_active_levels();
        auto& mask_cell = sim.repo().get_int_field("mask_cell");
        auto& velocity = sim.repo().get_field("velocity");

        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                auto epsilon_cell = mask_cell(lev).array(mfi);
                auto varr = velocity(lev).array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        if (epsilon_cell(i, j, k) == 0) {
                            varr(i, j, k, 0) = wdata.vel_bc[0];
                            varr(i, j, k, 1) = wdata.vel_bc[1];
                            varr(i, j, k, 2) = wdata.vel_bc[2];
                        }
                    });
            }
        }
    }
};

template <typename GeomTrait>
struct ComputeForceOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
{
    void operator()(typename GeomTrait::DataType&) {}
};

template <typename GeomTrait>
struct ProcessOutputsOp<
    GeomTrait,
    typename std::enable_if<
        std::is_base_of<BluffBodyType, GeomTrait>::value>::type>
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
        bluff_body::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        bluff_body::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), time.new_time());
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* BLUFF_BODY_OPS_H */
