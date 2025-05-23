#include "amr-wind/utilities/sampling/FreeSurfaceSampler.H"
#include "amr-wind/utilities/io_utils.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/index_operations.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

FreeSurfaceSampler::FreeSurfaceSampler(CFDSim& sim)
    : m_sim(sim), m_vof(sim.repo().get_field("vof"))
{}

FreeSurfaceSampler::~FreeSurfaceSampler() = default;

void FreeSurfaceSampler::initialize(const std::string& key)
{
    BL_PROFILE("amr-wind::FreeSurfaceSampler::initialize");

    {
        amrex::ParmParse pp(key);
        pp.getarr("plane_start", m_start);
        pp.getarr("plane_end", m_end);
        pp.getarr("plane_num_points", m_npts_dir);
        pp.query("search_direction", m_coorddir);
        pp.query("num_instances", m_ninst);
        pp.query("max_sample_points_per_cell", m_ncmax);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_start.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_end.size()) == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);
        check_bounds();

        switch (m_coorddir) {
        case 0: {
            m_gc0 = 1;
            m_gc1 = 2;
            break;
        }
        case 1: {
            m_gc0 = 0;
            m_gc1 = 2;
            break;
        }
        case 2: {
            m_gc0 = 0;
            m_gc1 = 1;
            break;
        }
        default: {
            amrex::Abort(
                "FreeSurfaceSampler: Invalid coordinate search direction "
                "encountered");
            break;
        }
        }
    }
    // Small number for floating-point comparisons
    constexpr amrex::Real eps = 1.0e-16;

    // Calculate total number of points
    m_npts = m_npts_dir[0] * m_npts_dir[1];

    // Turn parameters into 2D grid
    m_grid_locs.resize(m_npts);
    m_out.resize(static_cast<long>(m_npts) * m_ninst);

    // Get size of sample grid spacing
    amrex::Real dxs0 =
        (m_end[m_gc0] - m_start[m_gc0]) / amrex::max(m_npts_dir[0] - 1, 1);
    amrex::Real dxs1 =
        (m_end[m_gc1] - m_start[m_gc1]) / amrex::max(m_npts_dir[1] - 1, 1);

    // Store locations
    int idx = 0;
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            // Initialize output values to 0.0
            for (int ni = 0; ni < m_ninst; ++ni) {
                m_out[idx * m_ninst + ni] = m_start[m_coorddir];
            }
            // Grid direction 1
            m_grid_locs[idx][0] = m_start[m_gc0] + dxs0 * i;
            // Grid direction 2
            m_grid_locs[idx][1] = m_start[m_gc1] + dxs1 * j;

            ++idx;
        }
    }

    // Capture variables for device
    const amrex::Real s_gc0 = m_start[m_gc0];
    const amrex::Real s_gc1 = m_start[m_gc1];
    const int ntps0 = m_npts_dir[0];
    const int ntps1 = m_npts_dir[1];
    const int gc0 = m_gc0;
    const int gc1 = m_gc1;

    // Determine number of components necessary for working fields
    int ncomp = 0;
    const int finest_level = m_vof.repo().num_active_levels() - 1;
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), m_sim.mesh().refRatio(lev), 1,
                0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }
        // Get geometry information
        const auto& geom = m_sim.mesh().Geom(lev);
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
            geom.CellSizeArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
            geom.ProbLoArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> phi =
            geom.ProbHiArray();
        ncomp = amrex::max(
            ncomp,
            amrex::ReduceMax(
                level_mask, 0,
                [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Box const& bx,
                    amrex::Array4<int const> const& mask_arr) -> int {
                    int ns_fab = 0;
                    amrex::Loop(bx, [=, &ns_fab](int i, int j, int k) noexcept {
                        // Cell location
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
                        xm[0] = plo[0] + (i + 0.5) * dx[0];
                        xm[1] = plo[1] + (j + 0.5) * dx[1];
                        xm[2] = plo[2] + (k + 0.5) * dx[2];
                        int n0_f = 0;
                        int n0_a = 0;
                        int n1_f = 0;
                        int n1_a = 0;
                        // Get first and after sample indices for gc0
                        if (ntps0 == 1) {
                            n0_a = ((std::abs(phi[gc0] - s_gc0) < eps) ||
                                    (xm[gc0] - s_gc0 <= 0.5 * dx[gc0] &&
                                     s_gc0 - xm[gc0] < 0.5 * dx[gc0]))
                                       ? 1
                                       : 0;
                        } else {
                            n0_f = (int)std::ceil(
                                (xm[gc0] - 0.5 * dx[gc0] - s_gc0) / dxs0);
                            n0_a = (int)std::ceil(
                                (xm[gc0] + 0.5 * dx[gc0] - s_gc0) / dxs0);
                            // Edge case of phi
                            if (std::abs(xm[gc0] + 0.5 * dx[gc0] - phi[gc0]) <
                                    eps &&
                                std::abs(s_gc0 + n0_a * dxs0 - phi[gc0]) <
                                    eps) {
                                ++n0_a;
                            }
                            // Bounds
                            n0_a = amrex::min(ntps0, n0_a);
                            n0_f = amrex::max(amrex::min(0, n0_a), n0_f);
                            // Out of bounds indicates no sample point
                            if (n0_f >= ntps0 || n0_f < 0) {
                                n0_a = n0_f;
                            }
                        }
                        // Get first and after sample indices for gc1
                        if (ntps1 == 1) {
                            n1_a = (std::abs(phi[gc1] - s_gc1) < eps ||
                                    (xm[gc1] - s_gc1 <= 0.5 * dx[gc1] &&
                                     s_gc1 - xm[gc1] < 0.5 * dx[gc1]))
                                       ? 1
                                       : 0;
                        } else {
                            n1_f = (int)std::ceil(
                                (xm[gc1] - 0.5 * dx[gc1] - s_gc1) / dxs1);
                            n1_a = (int)std::ceil(
                                (xm[gc1] + 0.5 * dx[gc1] - s_gc1) / dxs1);
                            // Edge case of phi
                            if (std::abs(xm[gc1] + 0.5 * dx[gc1] - phi[gc1]) <
                                    eps &&
                                std::abs(s_gc1 + n1_a * dxs1 - phi[gc1]) <
                                    eps) {
                                ++n1_a;
                            }
                            // Bounds
                            n1_a = amrex::min(ntps1, n1_a);
                            n1_f = amrex::max(amrex::min(0, n1_a), n1_f);
                            // Out of bounds indicates no sample point
                            if (n1_f >= ntps1 || n1_f < 0) {
                                n1_a = n1_f;
                            }
                        }
                        // Get total number of possible samples in cell
                        int ns = (n0_a - n0_f) * (n1_a - n1_f);
                        ns_fab = amrex::max(ns_fab, mask_arr(i, j, k) * ns);
                    });
                    return ns_fab;
                }));
    }
    amrex::ParallelDescriptor::ReduceIntMax(ncomp);
    // Limit number of components
    ncomp = amrex::min(ncomp, m_ncmax);
    // Save number of components
    m_ncomp = ncomp;

    // Declare fields for search
    auto& floc =
        m_sim.repo().declare_field("sample_loc_" + m_label, 2 * ncomp, 0, 1);
    auto& fidx =
        m_sim.repo().declare_int_field("sample_idx_" + m_label, ncomp, 0, 1);

    // Store locations and indices in fields
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), m_sim.mesh().refRatio(lev), 1,
                0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }
        // Get geometry information
        const auto& geom = m_sim.mesh().Geom(lev);
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
            geom.CellSizeArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
            geom.ProbLoArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> phi =
            geom.ProbHiArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi(floc(lev)); mfi.isValid(); ++mfi) {
            auto loc_arr = floc(lev).array(mfi);
            auto idx_arr = fidx(lev).array(mfi);
            auto mask_arr = level_mask.const_array(mfi);
            const auto& vbx = mfi.validbox();
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Cell location
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
                    xm[0] = plo[0] + (i + 0.5) * dx[0];
                    xm[1] = plo[1] + (j + 0.5) * dx[1];
                    xm[2] = plo[2] + (k + 0.5) * dx[2];
                    int n0_f = 0;
                    int n0_a = 0;
                    int n1_f = 0;
                    int n1_a = 0;
                    // Get first and after sample indices for gc0
                    if (ntps0 == 1) {
                        n0_a = (std::abs(phi[gc0] - s_gc0) < eps ||
                                (xm[gc0] - s_gc0 <= 0.5 * dx[gc0] &&
                                 s_gc0 - xm[gc0] < 0.5 * dx[gc0]))
                                   ? 1
                                   : 0;
                    } else {
                        n0_f = (int)std::ceil(
                            (xm[gc0] - 0.5 * dx[gc0] - s_gc0) / dxs0);
                        n0_a = (int)std::ceil(
                            (xm[gc0] + 0.5 * dx[gc0] - s_gc0) / dxs0);
                        // Edge case of phi
                        if (std::abs(xm[gc0] + 0.5 * dx[gc0] - phi[gc0]) <
                                eps &&
                            std::abs(s_gc0 + n0_a * dxs0 - phi[gc0]) < eps) {
                            ++n0_a;
                        }
                        // Bounds
                        n0_a = amrex::min(ntps0, n0_a);
                        n0_f = amrex::max(amrex::min(0, n0_a), n0_f);
                        // Out of bounds indicates no sample point
                        if (n0_f >= ntps0 || n0_f < 0) {
                            n0_a = n0_f;
                        }
                    }
                    // Get first and after sample indices for gc1
                    if (ntps1 == 1) {
                        n1_a = (std::abs(phi[gc1] - s_gc1) < eps ||
                                (xm[gc1] - s_gc1 <= 0.5 * dx[gc1] &&
                                 s_gc1 - xm[gc1] < 0.5 * dx[gc1]))
                                   ? 1
                                   : 0;
                    } else {
                        n1_f = (int)std::ceil(
                            (xm[gc1] - 0.5 * dx[gc1] - s_gc1) / dxs1);
                        n1_a = (int)std::ceil(
                            (xm[gc1] + 0.5 * dx[gc1] - s_gc1) / dxs1);
                        // Edge case of phi
                        if (std::abs(xm[gc1] + 0.5 * dx[gc1] - phi[gc1]) <
                                eps &&
                            std::abs(s_gc1 + n1_a * dxs1 - phi[gc1]) < eps) {
                            ++n1_a;
                        }
                        // Bounds
                        n1_a = amrex::min(ntps1, n1_a);
                        n1_f = amrex::max(amrex::min(0, n1_a), n1_f);
                        // Out of bounds indicates no sample point
                        if (n1_f >= ntps1 || n1_f < 0) {
                            n1_a = n1_f;
                        }
                    }
                    // Loop through local sample locations
                    int ns = 0;
                    for (int n0 = n0_f; n0 < n0_a; ++n0) {
                        for (int n1 = n1_f; n1 < n1_a; ++n1) {
                            // Save index and location
                            idx_arr(i, j, k, ns) = n1 * ntps0 + n0;
                            loc_arr(i, j, k, 2 * ns) = s_gc0 + n0 * dxs0;
                            loc_arr(i, j, k, 2 * ns + 1) = s_gc1 + n1 * dxs1;
                            // Advance to next point
                            ++ns;
                            // if ns gets to max components, break
                            if (ns == ncomp) {
                                break;
                            }
                        }
                        if (ns == ncomp) {
                            break;
                        }
                    }
                    // Set remaining values to -1 to indicate no point
                    // or set all values to -1 if not in fine mesh
                    int nstart = (mask_arr(i, j, k) == 0) ? 0 : ns;
                    for (int n = nstart; n < ncomp; ++n) {
                        idx_arr(i, j, k, n) = -1;
                    }
                });
        }
    }
}
void FreeSurfaceSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    bool all_ok = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_start[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_start[d] = prob_lo[d] + 10 * bounds_tol;
        }
        if (m_start[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_start[d] = prob_hi[d] - 10 * bounds_tol;
        }
        if (m_end[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_end[d] = prob_lo[d] + 10 * bounds_tol;
        }
        if (m_end[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_end[d] = prob_hi[d] - 10 * bounds_tol;
        }
    }
    if (!all_ok) {
        amrex::Print()
            << "WARNING: FreeSurfaceSampler: Out of domain plane was "
               "truncated to match domain"
            << std::endl;
    }
}

void FreeSurfaceSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto domain = m_sim.mesh().Geom(lev).Domain();
    sampling_locations(sample_locs, domain);

    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void FreeSurfaceSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    int idx = 0;
    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    for (int j = 0; j < m_npts_dir[1]; ++j) {
        for (int i = 0; i < m_npts_dir[0]; ++i) {
            for (int ni = 0; ni < m_ninst; ++ni) {
                amrex::RealVect loc;
                loc[m_gc0] = m_grid_locs[idx][0];
                loc[m_gc1] = m_grid_locs[idx][1];
                loc[m_coorddir] = m_out[idx * m_ninst + ni];
                if (utils::contains(box, loc, plo, dxinv)) {
                    sample_locs.push_back(loc, idx * m_ninst + ni);
                }
            }
            ++idx;
        }
    }
}

bool FreeSurfaceSampler::update_sampling_locations()
{

    BL_PROFILE("amr-wind::FreeSurfaceSampler::update_sampling_locations");

    // Zero data in output array
    for (int n = 0; n < m_npts * m_ninst; n++) {
        m_out[n] = 0.0;
    }
    // Set up device vector of current outputs, initialize to plo
    const auto& plo0 = m_sim.mesh().Geom(0).ProbLoArray();
    amrex::Gpu::DeviceVector<amrex::Real> dout(m_npts, plo0[m_coorddir]);
    auto* dout_ptr = dout.data();
    // Set up device vector of last outputs, initialize to above phi0
    const auto& phi0 = m_sim.mesh().Geom(0).ProbHiArray();
    amrex::Gpu::DeviceVector<amrex::Real> dout_last(
        m_npts, phi0[m_coorddir] + 1.0);
    auto* dlst_ptr = dout_last.data();

    // Get working fields
    auto& fidx = m_sim.repo().get_int_field("sample_idx_" + m_label);
    auto& floc = m_sim.repo().get_field("sample_loc_" + m_label);

    const int finest_level = m_vof.repo().num_active_levels() - 1;

    // Capture integers for device
    const int dir = m_coorddir;
    const int gc0 = m_gc0;
    const int gc1 = m_gc1;
    const int ncomp = m_ncomp;

    // Loop instances
    for (int ni = 0; ni < m_ninst; ++ni) {
        for (int lev = 0; lev <= finest_level; lev++) {
            // Level mask info is built into idx info

            // Get geometry information
            const auto& geom = m_sim.mesh().Geom(lev);
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
                geom.CellSizeArray();
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxi =
                geom.InvCellSizeArray();
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
                geom.ProbLoArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
            for (amrex::MFIter mfi(floc(lev)); mfi.isValid(); ++mfi) {
                auto loc_arr = floc(lev).const_array(mfi);
                auto idx_arr = fidx(lev).const_array(mfi);
                auto vof_arr = m_vof(lev).const_array(mfi);
                const auto& vbx = mfi.validbox();
                amrex::ParallelFor(
                    vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Cell location
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
                        xm[0] = plo[0] + (i + 0.5) * dx[0];
                        xm[1] = plo[1] + (j + 0.5) * dx[1];
                        xm[2] = plo[2] + (k + 0.5) * dx[2];
                        // Loop number of components
                        for (int n = 0; n < ncomp; ++n) {
                            // Get index of current component and cell
                            const int idx = idx_arr(i, j, k, n);
                            // Proceed if there is sample point at this i,j,k,n
                            // and that cell height is below previous instance
                            if (idx >= 0 && dlst_ptr[amrex::max(0, idx)] >
                                                xm[dir] + 0.5 * dx[dir]) {
                                // Get sample locations from field loc
                                amrex::Real loc0 = loc_arr(i, j, k, 2 * n);
                                amrex::Real loc1 = loc_arr(i, j, k, 2 * n + 1);

                                // Indices and slope variables
                                // Get modified indices for checking up
                                // and down and orient normal in search
                                // direction
                                amrex::Real mx = (dir == 0) ? 1.0 : 0.0;
                                amrex::Real my = (dir == 1) ? 1.0 : 0.0;
                                amrex::Real mz = (dir == 2) ? 1.0 : 0.0;
                                amrex::Real alpha = 1.0;
                                // If cell is full of single phase
                                // (accounts for when interface is at
                                // intersection of cells but lower one
                                // is not single-phase)
                                bool calc_flag = false;
                                if ((ni % 2 == 0 &&
                                     vof_arr(i, j, k) >= 1.0 - 1e-12) ||
                                    (ni % 2 != 0 &&
                                     vof_arr(i, j, k) <= 1e-12)) {
                                    // put bdy at top
                                    alpha = 1.0;
                                    if (ni % 2 != 0) {
                                        mx *= -1.0;
                                        my *= -1.0;
                                        mz *= -1.0;
                                        alpha *= -1.0;
                                    }
                                    calc_flag = true;
                                }
                                // Multiphase cell case
                                if (vof_arr(i, j, k) < (1.0 - 1e-12) &&
                                    vof_arr(i, j, k) > 1e-12) {
                                    amrex::IntVect iv_up{i, j, k},
                                        iv_down{i, j, k};
                                    iv_up[dir] += 1;
                                    iv_down[dir] -= 1;
                                    const bool intersect_above =
                                        (vof_arr(i, j, k) - 0.5) *
                                            (vof_arr(iv_up) - 0.5) <=
                                        0.;
                                    const bool intersect_below =
                                        (vof_arr(i, j, k) - 0.5) *
                                            (vof_arr(iv_down) - 0.5) <=
                                        0.;
                                    const bool closer_than_above =
                                        std::abs(vof_arr(i, j, k) - 0.5) <
                                        std::abs(vof_arr(iv_up) - 0.5);
                                    const bool closer_than_below =
                                        std::abs(vof_arr(i, j, k) - 0.5) <=
                                        std::abs(vof_arr(iv_down) - 0.5);
                                    calc_flag =
                                        (intersect_above &&
                                         closer_than_above) ||
                                        (intersect_below && closer_than_below);
                                    // Get interface reconstruction
                                    if (calc_flag) {
                                        multiphase::fit_plane(
                                            i, j, k, vof_arr, mx, my, mz,
                                            alpha);
                                    }
                                }

                                if (calc_flag) {
                                    // Initialize height measurement
                                    amrex::Real ht = plo[dir];
                                    // Reassign slope coefficients
                                    const amrex::Real mdr =
                                        (dir == 0) ? mx
                                                   : ((dir == 1) ? my : mz);
                                    const amrex::Real mg1 =
                                        (dir == 0) ? my : mx;
                                    const amrex::Real mg2 =
                                        (dir == 2) ? my : mz;
                                    // Get height of interface
                                    if (mdr == 0) {
                                        // If slope is undefined in z,
                                        // use middle of cell
                                        ht = xm[dir];
                                    } else {
                                        // Intersect 2D point with plane
                                        ht = (xm[dir] - 0.5 * dx[dir]) +
                                             (alpha -
                                              mg1 * dxi[gc0] *
                                                  (loc0 -
                                                   (xm[gc0] - 0.5 * dx[gc0])) -
                                              mg2 * dxi[gc1] *
                                                  (loc1 -
                                                   (xm[gc1] - 0.5 * dx[gc1]))) /
                                                 (mdr * dxi[dir]);
                                    }
                                    // If interface is below lower
                                    // bound, continue to look
                                    if (ht < xm[dir] - 0.5 * dx[dir]) {
                                        ht = plo[dir];
                                    }
                                    // If interface is above upper
                                    // bound, limit it
                                    if (ht > xm[dir] +
                                                 0.5 * dx[dir] * (1.0 + 1e-8)) {
                                        ht = xm[dir] + 0.5 * dx[dir];
                                    }
                                    // Save interface location by atomic max
                                    amrex::Gpu::Atomic::Max(&dout_ptr[idx], ht);
                                }
                            }
                        }
                    });
            }
        }

        // Copy information back from device
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, dout.begin(), dout.end(),
            &m_out[static_cast<long>(ni) * m_npts]);
        // Make consistent across parallelization
        for (int n = 0; n < m_npts; n++) {
            amrex::ParallelDescriptor::ReduceRealMax(m_out[ni * m_npts + n]);
        }
        // Copy last m_out to device vector of results of last instance
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, &m_out[static_cast<long>(ni) * m_npts],
            &m_out[(ni + 1) * m_npts - 1] + 1, dout_last.begin());
        // Reset current output device vector
        const auto coorddir = m_coorddir;
        amrex::ParallelFor(m_npts, [=] AMREX_GPU_DEVICE(int n) {
            dout_ptr[n] = plo0[coorddir];
        });
    }

    return true;
}

void FreeSurfaceSampler::post_regrid_actions()
{
    BL_PROFILE("amr-wind::FreeSurfaceSampler::post_regrid_actions");
    // Small number for floating-point comparisons
    constexpr amrex::Real eps = 1.0e-16;
    // Get working fields
    auto& fidx = m_sim.repo().get_int_field("sample_idx_" + m_label);
    auto& floc = m_sim.repo().get_field("sample_loc_" + m_label);
    // Provide variables from class to device
    const int ncomp = m_ncomp;
    const int ntps0 = m_npts_dir[0];
    const int ntps1 = m_npts_dir[1];
    const int gc0 = m_gc0;
    const int gc1 = m_gc1;
    const amrex::Real s_gc0 = m_start[m_gc0];
    const amrex::Real s_gc1 = m_start[m_gc1];
    const amrex::Real dxs0 =
        (m_end[m_gc0] - m_start[m_gc0]) / amrex::max(m_npts_dir[0] - 1, 1);
    const amrex::Real dxs1 =
        (m_end[m_gc1] - m_start[m_gc1]) / amrex::max(m_npts_dir[1] - 1, 1);

    // Store locations and indices in fields
    const int finest_level = m_vof.repo().num_active_levels() - 1;
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), m_sim.mesh().refRatio(lev), 1,
                0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }
        // Get geometry information
        const auto& geom = m_sim.mesh().Geom(lev);
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
            geom.CellSizeArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo =
            geom.ProbLoArray();
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> phi =
            geom.ProbHiArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi(floc(lev)); mfi.isValid(); ++mfi) {
            auto loc_arr = floc(lev).array(mfi);
            auto idx_arr = fidx(lev).array(mfi);
            auto mask_arr = level_mask.const_array(mfi);
            const auto& vbx = mfi.validbox();
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Cell location
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xm;
                    xm[0] = plo[0] + (i + 0.5) * dx[0];
                    xm[1] = plo[1] + (j + 0.5) * dx[1];
                    xm[2] = plo[2] + (k + 0.5) * dx[2];
                    int n0_f = 0;
                    int n0_a = 0;
                    int n1_f = 0;
                    int n1_a = 0;
                    // Get first and after sample indices for gc0
                    if (ntps0 == 1) {
                        n0_a = (std::abs(phi[gc0] - s_gc0) < eps ||
                                (xm[gc0] - s_gc0 <= 0.5 * dx[gc0] &&
                                 s_gc0 - xm[gc0] < 0.5 * dx[gc0]))
                                   ? 1
                                   : 0;
                    } else {
                        n0_f = (int)std::ceil(
                            (xm[gc0] - 0.5 * dx[gc0] - s_gc0) / dxs0);
                        n0_a = (int)std::ceil(
                            (xm[gc0] + 0.5 * dx[gc0] - s_gc0) / dxs0);
                        // Edge case of phi
                        if (std::abs(xm[gc0] + 0.5 * dx[gc0] - phi[gc0]) <
                                eps &&
                            std::abs(s_gc0 + n0_a * dxs0 - phi[gc0]) < eps) {
                            ++n0_a;
                        }
                        // Bounds
                        n0_a = amrex::min(ntps0, n0_a);
                        n0_f = amrex::max(amrex::min(0, n0_a), n0_f);
                        // Out of bounds indicates no sample point
                        if (n0_f >= ntps0 || n0_f < 0) {
                            n0_a = n0_f;
                        }
                    }
                    // Get first and after sample indices for gc1
                    if (ntps1 == 1) {
                        n1_a = (std::abs(phi[gc1] - s_gc1) < eps ||
                                (xm[gc1] - s_gc1 <= 0.5 * dx[gc1] &&
                                 s_gc1 - xm[gc1] < 0.5 * dx[gc1]))
                                   ? 1
                                   : 0;
                    } else {
                        n1_f = (int)std::ceil(
                            (xm[gc1] - 0.5 * dx[gc1] - s_gc1) / dxs1);
                        n1_a = (int)std::ceil(
                            (xm[gc1] + 0.5 * dx[gc1] - s_gc1) / dxs1);
                        // Edge case of phi
                        if (std::abs(xm[gc1] + 0.5 * dx[gc1] - phi[gc1]) <
                                eps &&
                            std::abs(s_gc1 + n1_a * dxs1 - phi[gc1]) < eps) {
                            ++n1_a;
                        }
                        // Bounds
                        n1_a = amrex::min(ntps1, n1_a);
                        n1_f = amrex::max(amrex::min(0, n1_a), n1_f);
                        // Out of bounds indicates no sample point
                        if (n1_f >= ntps1 || n1_f < 0) {
                            n1_a = n1_f;
                        }
                    }
                    // Loop through local sample locations
                    int ns = 0;
                    for (int n0 = n0_f; n0 < n0_a; ++n0) {
                        for (int n1 = n1_f; n1 < n1_a; ++n1) {
                            // Save index and location
                            idx_arr(i, j, k, ns) = n1 * ntps0 + n0;
                            loc_arr(i, j, k, 2 * ns) = s_gc0 + n0 * dxs0;
                            loc_arr(i, j, k, 2 * ns + 1) = s_gc1 + n1 * dxs1;
                            // Advance to next point
                            ++ns;
                            // if ns gets to max components, break
                            if (ns == ncomp) {
                                break;
                            }
                        }
                        if (ns == ncomp) {
                            break;
                        }
                    }
                    // Set remaining values to -1 to indicate no point
                    // or set all values to -1 if not in fine mesh
                    int nstart = (mask_arr(i, j, k) == 0) ? 0 : ns;
                    for (int n = nstart; n < ncomp; ++n) {
                        idx_arr(i, j, k, n) = -1;
                    }
                });
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void FreeSurfaceSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{

    // Metadata related to the 2D grid used to sample
    const std::vector<int> ijk{m_npts_dir[0], m_npts_dir[1], m_ninst};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("plane_start", m_start);
    grp.put_attr("plane_end", m_end);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void FreeSurfaceSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void FreeSurfaceSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    // Write the coordinates every time
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    SampleLocType sample_locs;
    sampling_locations(sample_locs);
    auto xyz = grp.var("points");
    count[1] = num_output_points();
    const auto& locs = sample_locs.locations();
    xyz.put(locs[0].begin(), start, count);
}
#else
void FreeSurfaceSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void FreeSurfaceSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void FreeSurfaceSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
