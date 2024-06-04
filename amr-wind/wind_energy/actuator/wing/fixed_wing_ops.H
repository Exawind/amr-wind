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

        amrex::Vector<amrex::Real> fcf, dcf;
        // If certain components of Gaussian should be ignored
        pp.queryarr("active_gaussian_dirs", dcf);
        // If only a force in fx direction should be applied
        pp.queryarr("active_force_dirs", fcf);
        // Populate vectors in data
        if (!dcf.empty()) {
            if (dcf.size() == AMREX_SPACEDIM) {
                wdata.dcoord_flags.x() = dcf[0];
                wdata.dcoord_flags.y() = dcf[1];
                wdata.dcoord_flags.z() = dcf[2];
            } else {
                amrex::Abort(
                    "fixed_wing: wrong number of components specified for "
                    "active_gaussian_dirs.\n");
            }
        }
        if (!fcf.empty()) {
            if (fcf.size() == AMREX_SPACEDIM) {
                wdata.fcoord_flags.x() = fcf[0];
                wdata.fcoord_flags.y() = fcf[1];
                wdata.fcoord_flags.z() = fcf[2];
            } else {
                amrex::Abort(
                    "fixed_wing: wrong number of components specified for "
                    "active_force_dirs.\n");
            }
        }

        // If force magnitude should be computed with |U|=1
        pp.query("uinf_unity", wdata.uinf_unity);

        // Initialize tables for pitch actuation (in degrees)
        pp.query("pitch_timetable", wdata.pitch_timetable_file);
        if (!wdata.pitch_timetable_file.empty()) {
            std::ifstream ifh(wdata.pitch_timetable_file, std::ios::in);
            if (!ifh.good()) {
                amrex::Abort(
                    "Cannot find input file: " + wdata.pitch_timetable_file);
            }
            amrex::Real data_time;
            amrex::Real data_pitch_deg;
            ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            while (ifh >> data_time) {
                ifh >> data_pitch_deg;
                wdata.time_table.push_back(data_time);
                wdata.pitch_table.push_back(data_pitch_deg);
            }
        }

        pp.query("motion_type", wdata.motion_type);
        {
            // do nothing for default, "none"
            if (amrex::toLower(wdata.motion_type) == "linear") {
                // get velocity
                pp.get("velocity", wdata.vel_tr);
            } else if (amrex::toLower(wdata.motion_type) == "sine") {
                // get parameters
                pp.get("sine_period", wdata.s_period);
                pp.get("sine_vector", wdata.s_vector);
            } else if (amrex::toLower(wdata.motion_type) != "none") {
                amrex::Abort(
                    "Invalid FixedWingLine motion type. Valid options are "
                    "none, linear, and sine.");
            }
        }
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
        // Set up bounding box depending on motion type
        if (amrex::toLower(wdata.motion_type) == "none") {
            // clang-format off
            info.bound_box = amrex::RealBox(
                amrex::min(p1.x(), p2.x()) - search_radius,
                amrex::min(p1.y(), p2.y()) - search_radius,
                amrex::min(p1.z(), p2.z()) - search_radius,
                amrex::max(p1.x(), p2.x()) + search_radius,
                amrex::max(p1.y(), p2.y()) + search_radius,
                amrex::max(p1.z(), p2.z()) + search_radius);
            // clang-format on
        } else if (amrex::toLower(wdata.motion_type) == "linear") {
            // Extend bounding box in case of velocity
            constexpr amrex::Real tiny = 1e-8;
            // If velocity is present, assume wing can travel entire domain
            const bool up_ = wdata.vel_tr.x() > tiny;
            const amrex::Real up_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[0];
            const bool um_ = wdata.vel_tr.x() < -tiny;
            const amrex::Real um_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[0];
            const bool vp_ = wdata.vel_tr.y() > tiny;
            const amrex::Real vp_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[1];
            const bool vm_ = wdata.vel_tr.y() < -tiny;
            const amrex::Real vm_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[1];
            const bool wp_ = wdata.vel_tr.z() > tiny;
            const amrex::Real wp_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[2];
            const bool wm_ = wdata.vel_tr.z() < -tiny;
            const amrex::Real wm_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[2];
            info.bound_box = amrex::RealBox(
                um_ ? um_ext : amrex::min(p1.x(), p2.x()) - search_radius,
                vm_ ? vm_ext : amrex::min(p1.y(), p2.y()) - search_radius,
                wm_ ? wm_ext : amrex::min(p1.z(), p2.z()) - search_radius,
                up_ ? up_ext : amrex::max(p1.x(), p2.x()) + search_radius,
                vp_ ? vp_ext : amrex::max(p1.y(), p2.y()) + search_radius,
                wp_ ? wp_ext : amrex::max(p1.z(), p2.z()) + search_radius);
        } else if (amrex::toLower(wdata.motion_type) == "sine") {
            // clang-format off
            info.bound_box = amrex::RealBox(
                amrex::min(p1.x(), p2.x()) - wdata.s_vector.x() - search_radius,
                amrex::min(p1.y(), p2.y()) - wdata.s_vector.y() - search_radius,
                amrex::min(p1.z(), p2.z()) - wdata.s_vector.z() - search_radius,
                amrex::max(p1.x(), p2.x()) + wdata.s_vector.x() + search_radius,
                amrex::max(p1.y(), p2.y()) + wdata.s_vector.y() + search_radius,
                amrex::max(p1.z(), p2.z()) + wdata.s_vector.z() + search_radius);
            // clang-format on
        }
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
                meta.epsilon_chord.x(), meta.epsilon_chord.y(),
                meta.epsilon_chord.z()};
            // clang-format on
            for (int i = 0; i < npts; ++i) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    const auto eps = meta.epsilon_chord[n] * meta.chord[i];
                    grid.epsilon[i][n] = amrex::max(meta.eps_inp[n], eps);
                }
            }
            // Copy Gaussian flags to grid struct
            grid.dcoord_flags = meta.dcoord_flags;
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
