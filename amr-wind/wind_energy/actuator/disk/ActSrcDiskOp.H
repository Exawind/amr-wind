#ifndef ACTSRCDISKOP_H_
#define ACTSRCDISKOP_H_
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/wind_energy/actuator/disk/disk_types.H"

namespace amr_wind {
namespace actuator {
namespace ops {

template <typename ActTrait>
class ActSrcOp<
    ActTrait,
    ActSrcDisk,
    typename std::enable_if<std::is_base_of<DiskType, ActTrait>::value>::type>
{
private:
    typename ActTrait::DataType& m_data;
    Field& m_act_src;

    DeviceVecList m_pos;
    DeviceVecList m_force;

    void copy_to_device();

public:
    explicit ActSrcOp(typename ActTrait::DataType& data)
        : m_data(data)
        , m_act_src(m_data.sim().repo().get_field("actuator_src_term"))
    {}

    void initialize();

    void setup_op() { copy_to_device(); }

    void operator()(
        const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);
};

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    typename std::enable_if<std::is_base_of<DiskType, ActTrait>::value>::type>::
    initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_force.resize(grid.force.size());
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    typename std::enable_if<std::is_base_of<DiskType, ActTrait>::value>::type>::
    copy_to_device()
{
    const auto& grid = m_data.grid();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.force.begin(), grid.force.end(),
        m_force.begin());
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    typename std::enable_if<std::is_base_of<DiskType, ActTrait>::value>::type>::
operator()(const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    const std::string fname = ActTrait::identifier() + ActSrcDisk::identifier();
    BL_PROFILE("amr-wind::ActSrcOp<" + fname + ">");

    const auto& bx = mfi.tilebox();
    const auto& sarr = m_act_src(lev).array(mfi);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    auto& data = m_data.meta();

    const amrex::Real dR = data.dr;
    const amrex::Real epsilon = data.epsilon;
    const vs::Vector m_normal(data.normal_vec);
    const vs::Vector m_origin(data.center);
    const auto* pos = m_pos.data();
    const auto* force = m_force.data();
    const int npts = data.num_force_pts;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const vs::Vector cc{
            problo[0] + (i + 0.5) * dx[0],
            problo[1] + (j + 0.5) * dx[1],
            problo[2] + (k + 0.5) * dx[2],
        };

        amrex::Real src_force[AMREX_SPACEDIM]{0.0, 0.0, 0.0};
        for (int ip = 0; ip < npts; ++ip) {
            const auto R =
                utils::delta_pnts_cyl(m_origin, m_normal, m_origin, pos[ip])
                    .x();
            const auto dist_on_disk =
                utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
            const auto& pforce = force[ip];

            const amrex::Real weight_R =
                utils::linear_basis_1d(dist_on_disk.x(), dR);
            const amrex::Real weight_T =
                1.0 / (::amr_wind::utils::two_pi() * R);
            const amrex::Real weight_N =
                utils::gaussian1d(dist_on_disk.z(), epsilon);
            const auto projection_weight = weight_R * weight_T * weight_N;

            src_force[0] += projection_weight * pforce.x();
            src_force[1] += projection_weight * pforce.y();
            src_force[2] += projection_weight * pforce.z();
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}
} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* ACTSRCDISKOP_H_ */
