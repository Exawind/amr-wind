#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/geometry/Sphere.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Sphere>
{
    void operator()(Sphere::DataType& data, const utils::IBParser& pp)
    {
        auto& wdata = data.meta();
        pp.get("num_points", wdata.num_pts);
        // clang-format on
    }
};

template <>
struct InitDataOp<Sphere>
{
    void operator()(Sphere::DataType&)
    {
        /*
            wing::init_data_structures(data.meta(), data.grid());

            auto& meta = data.meta();
            {
                const int npts = data.meta().num_pts;
                auto& grid = data.grid();
                // Process chord information
                const auto wlen = vs::mag(grid.pos.back() - grid.pos.front());
                RealList wx(npts);
                for (int i = 0; i < npts; ++i) {
                    wx[i] = vs::mag(grid.pos[i]) / wlen;
                }
                meta.chord.resize(npts);
                ::amr_wind::interp::linear_monotonic(
                    meta.span_locs, meta.chord_inp, wx, meta.chord);
                // clang-format off
                meta.epsilon_chord = {
                    meta.epsilon_chord.x(), meta.epsilon_chord.z(),
                    meta.epsilon_chord.y()};
                // clang-format on
                for (int i = 0; i < npts; ++i) {
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        const auto eps = meta.epsilon_chord[n] * meta.chord[i];
                        grid.epsilon[i][n] = amrex::max(meta.eps_inp[n], eps);
                    }
                }
            }
        }
        */
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* SPHERE_OPS_H */
