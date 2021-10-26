#ifndef ACTSRCDISKOP_H_
#define ACTSRCDISKOP_H_
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/wind_energy/actuator/disk/disk_types.H"

namespace amr_wind {
namespace actuator {
namespace ops {

template <typename ActTrait, typename SpreadingType>
class ActSrcOp<
    ActTrait,
    ActSrcDiskBase<SpreadingType>,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>
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

    template <typename U = SpreadingType>
    std::enable_if_t<std::is_same<U, LinearBasis>::value> operator()(
        const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
    {
        const std::string fname =
            ActTrait::identifier() + ActSrcDisk::identifier();
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

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const vs::Vector cc{
                    problo[0] + (i + 0.5) * dx[0],
                    problo[1] + (j + 0.5) * dx[1],
                    problo[2] + (k + 0.5) * dx[2],
                };

                amrex::Real src_force[AMREX_SPACEDIM]{0.0, 0.0, 0.0};
                for (int ip = 0; ip < npts; ++ip) {
                    const auto R = utils::delta_pnts_cyl(
                                       m_origin, m_normal, m_origin, pos[ip])
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
                    const auto projection_weight =
                        weight_R * weight_T * weight_N;

                    src_force[0] += projection_weight * pforce.x();
                    src_force[1] += projection_weight * pforce.y();
                    src_force[2] += projection_weight * pforce.z();
                }

                sarr(i, j, k, 0) += src_force[0];
                sarr(i, j, k, 1) += src_force[1];
                sarr(i, j, k, 2) += src_force[2];
            });
    }

    template <typename U = SpreadingType>
    std::enable_if_t<std::is_same<U, UniformGaussian>::value> operator()(
        const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
    {
        const std::string fname =
            ActTrait::identifier() + ActSrcDisk::identifier();
        BL_PROFILE("amr-wind::ActSrcOp<" + fname + ">");

        const auto& bx = mfi.tilebox();
        const auto& sarr = m_act_src(lev).array(mfi);
        const auto& problo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();

        auto& data = m_data.meta();

        const amrex::Real dR = data.dr;
        const vs::Vector epsilon = vs::Vector::one() * data.epsilon;
        const vs::Vector m_normal(data.normal_vec);
        const vs::Vector m_origin(data.center);
        const auto* pos = m_pos.data();
        const auto* force = m_force.data();
        const int npts = data.num_force_pts;
        const int nForceTheta = data.num_force_theta_pts;
        const auto dTheta = ::amr_wind::utils::two_pi() / nForceTheta;

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const vs::Vector cc{
                    problo[0] + (i + 0.5) * dx[0],
                    problo[1] + (j + 0.5) * dx[1],
                    problo[2] + (k + 0.5) * dx[2],
                };

                amrex::Real src_force[AMREX_SPACEDIM]{0.0, 0.0, 0.0};
                for (int ip = 0; ip < npts; ++ip) {
                    const auto dist_ref_pnt =
                        utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
                    const auto& pforce = force[ip];

                    for (int it = 0; it < nForceTheta; ++it) {
                        // loop over num points in theta and subtract theta
                        // change to get force point location
                        const auto cyl_dist =
                            dist_ref_pnt - vs::Vector(0, it * dTheta, 0);
                        // need to convert from cylindrical coordinates to
                        // cartesian distance is between points so we can use
                        // origin as second point
                        const vs::Vector cart_dist(
                            cyl_dist.x() * std::cos(cyl_dist.y()),
                            cyl_dist.x() * std::sin(cyl_dist.y()),
                            cyl_dist.z());

                        const auto projection_weight =
                            utils::gaussian3d(cart_dist, epsilon);

                        src_force[0] += projection_weight * pforce.x();
                        src_force[1] += projection_weight * pforce.y();
                        src_force[2] += projection_weight * pforce.z();
                    }
                }

                sarr(i, j, k, 0) += src_force[0];
                sarr(i, j, k, 1) += src_force[1];
                sarr(i, j, k, 2) += src_force[2];
            });
    }
};

template <typename ActTrait, typename SpreadingType>
void ActSrcOp<
    ActTrait,
    ActSrcDiskBase<SpreadingType>,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_force.resize(grid.force.size());
}

template <typename ActTrait, typename SpreadingType>
void ActSrcOp<
    ActTrait,
    ActSrcDiskBase<SpreadingType>,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>::
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

/*template <typename ActTrait, typename SpreadingType>
std::enable_if_t<std::is_same<SpreadingType, LinearBasis>::value> ActSrcOp<
    ActTrait,
    ActSrcDiskBase<SpreadingType>,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>::
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

template <typename ActTrait, typename SpreadingType>
std::enable_if_t<std::is_same<SpreadingType, UniformGaussian>::value> ActSrcOp<
    ActTrait,
    ActSrcDiskBase<SpreadingType>,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>::
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
    const int nForceTheta = data.num_force_theta_pts;
    const auto dTheta = ::amr_wind::tools::two_pi() / nForceTheta;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const vs::Vector cc{
            problo[0] + (i + 0.5) * dx[0],
            problo[1] + (j + 0.5) * dx[1],
            problo[2] + (k + 0.5) * dx[2],
        };

        amrex::Real src_force[AMREX_SPACEDIM]{0.0, 0.0, 0.0};
        for (int ip = 0; ip < npts; ++ip) {
            const auto dist_ref_pnt =
                utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
            const auto& pforce = force[ip];

            for (int it = 0; it < nForceTheta; ++it) {
                // loop over num points in theta and subtract theta
                // change to get force point location
                const auto cyl_dist =
                    dist_ref_pnt - vs::Vector(0, it * dTheta, 0);
                // need to convert from cylindrical coordinates to
                // cartesian distance is between points so we can use
                // origin as second point
                const vs::Vector cart_dist(
                    cyl_dist.x() * std::cos(cyl_dist.y()),
                    cyl_dist.x() * std::sin(cyl_dist.y()), cyl_dist.z());
                const double dist = vs::mag(cart_dist);

                const auto projection_weight = utils::gaussian3d(dist, epsilon);

                src_force[0] += projection_weight * pforce.x();
                src_force[1] += projection_weight * pforce.y();
                src_force[2] += projection_weight * pforce.z();
            }
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}*/
} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* ACTSRCDISKOP_H_ */
