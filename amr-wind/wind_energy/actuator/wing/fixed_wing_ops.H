#ifndef FIXED_WING_OPS_H
#define FIXED_WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/FixedWing.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::actuator::ops {

ComponentView make_component_view(FixedWing::DataType& data)
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

template <>
struct ReadInputsOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.query("epsilon", wdata.eps_inp);
        pp.query("epsilon_chord", wdata.epsilon_chord);
        pp.get("pitch", wdata.pitch);
        pp.get("airfoil_table", wdata.airfoil_file);
        pp.query("airfoil_type", wdata.airfoil_type);
        pp.queryarr("span_locs", wdata.span_locs);
        pp.queryarr("chord", wdata.chord_inp);
        bool use_fllc = false;
        pp.query("fllc", use_fllc);
        if (use_fllc) {
            wdata.fllc = std::make_unique<FLLCData>();
            FLLCParse(pp, *(wdata.fllc));
        }

        if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
            amrex::Abort(
                "Actuator fixed wing requires specification of one or both "
                "of 'epsilon' or 'epsilon_chord'");
        }

        AMREX_ALWAYS_ASSERT(wdata.span_locs.size() == wdata.chord_inp.size());

        amrex::Real max_chord =
            *std::max_element(wdata.chord_inp.begin(), wdata.chord_inp.end());
        amrex::Real max_eps =
            *std::max_element(wdata.eps_inp.begin(), wdata.eps_inp.end());
        amrex::Real max_epsc = *std::max_element(
            wdata.epsilon_chord.begin(), wdata.epsilon_chord.end());
        amrex::Real search_radius =
            amrex::max(max_eps, max_epsc) * max_chord * 3.0;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // clang-format off
        info.bound_box = amrex::RealBox(
            amrex::min(p1.x(), p2.x()) - search_radius,
            amrex::min(p1.y(), p2.y()) - search_radius,
            amrex::min(p1.z(), p2.z()) - search_radius,
            amrex::max(p1.x(), p2.x()) + search_radius,
            amrex::max(p1.y(), p2.y()) + search_radius,
            amrex::max(p1.z(), p2.z()) + search_radius
        );
        // clang-format on
    }
};

template <>
struct InitDataOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data)
    {
        wing::init_data_structures(data.meta(), data.grid());

        auto& meta = data.meta();
        {
            const int npts = data.meta().num_pts;
            auto& grid = data.grid();
            // Process chord information
            const auto wlen = vs::mag(grid.pos.back() - grid.pos.front());
            RealList wx(npts);
            for (int i = 0; i < npts; ++i) {
                wx[i] = vs::mag(grid.pos[i] - grid.pos[0]) / wlen;
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

        meta.aflookup =
            AirfoilLoader::load_airfoil(meta.airfoil_file, meta.airfoil_type);
    }
};

template <>
inline const AirfoilTable& airfoil_lookup<FixedWing>(FixedWing::DataType& data)
{
    return *data.meta().aflookup;
}

} // namespace amr_wind::actuator::ops

#endif /* FIXED_WING_OPS_H */
