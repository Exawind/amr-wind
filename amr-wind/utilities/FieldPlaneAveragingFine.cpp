#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include <algorithm>

namespace amr_wind {

template <typename FType>
FPlaneAveragingFine<FType>::FPlaneAveragingFine(
    const FType& field_in, const amr_wind::SimTime& time, int axis_in)
    : m_field(field_in), m_time(time), m_axis(axis_in)
{
    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);
    auto geom = m_field.repo().mesh().Geom();

    // beginning and end of line, for now assuming line is the length of the
    // entire domain
    m_xlo = geom[0].ProbLo(m_axis);
    m_xhi = geom[0].ProbHi(m_axis);

    const int finestLevel = m_field.repo().mesh().maxLevel();
    const int dom_lo = geom[finestLevel].Domain().smallEnd()[m_axis];
    const int dom_hi = geom[finestLevel].Domain().bigEnd()[m_axis];

    AMREX_ALWAYS_ASSERT(dom_lo == 0);
    int dom_hi2 = geom[0].Domain().bigEnd()[m_axis] + 1;
    for (int i = 0; i < finestLevel; ++i) {
        dom_hi2 *= 2;
    }
    AMREX_ALWAYS_ASSERT(dom_hi + 1 == dom_hi2);

    // FIXME: make an input maybe?
    m_ncell_line = dom_hi - dom_lo + 1;

    m_dx = (m_xhi - m_xlo) / static_cast<amrex::Real>(m_ncell_line);

    m_ncomp = m_field.num_comp();

    m_line_average.resize(static_cast<size_t>(m_ncell_line) * m_ncomp, 0.0);
    m_line_xcentroid.resize(m_ncell_line);

    for (int i = 0; i < m_ncell_line; ++i) {
        m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dx;
    }
}

template <typename FType>
void FPlaneAveragingFine<FType>::convert_x_to_ind(
    amrex::Real x, int& ind, amrex::Real& c) const
{
    c = 0.0;
    ind = 0;

    if (x > m_xlo + 0.5 * m_dx) {
        ind = static_cast<int>(floor((x - m_xlo) / m_dx - 0.5));
        const amrex::Real x1 = m_xlo + (ind + 0.5) * m_dx;
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);
}

template <typename FType>
void FPlaneAveragingFine<FType>::output_line_average_ascii(
    const std::string& filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::FPlaneAveragingFine::output_line_average_ascii");

    if (step != m_last_updated_index) {
        operator()();
    }

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    std::ofstream outfile;
    outfile.precision(m_precision);

    if (step == 1) {
        // make new file
        outfile.open(filename.c_str(), std::ios_base::out);
        outfile << "#ncell,ncomp" << std::endl;
        outfile << m_ncell_line << ", " << m_ncomp + 3 << std::endl;
        outfile << "#step,time,z";
        for (int i = 0; i < m_ncomp; ++i) {
            outfile << ",<" + m_field.name() + std::to_string(i) + ">";
        }
        outfile << std::endl;
    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    for (int i = 0; i < m_ncell_line; ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_line_xcentroid[i];
        for (int n = 0; n < m_ncomp; ++n) {
            outfile << ", " << std::scientific
                    << m_line_average[m_ncomp * i + n];
        }
        outfile << std::endl;
    }
}

template <typename FType>
void FPlaneAveragingFine<FType>::output_line_average_ascii(
    int step, amrex::Real time)
{
    const std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

template <typename FType>
amrex::Real FPlaneAveragingFine<FType>::line_average_interpolated(
    amrex::Real x, int comp) const
{

    BL_PROFILE("amr-wind::PlaneAveragingFine::line_average_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_average[m_ncomp * ind + comp] * (1.0 - c) +
           m_line_average[m_ncomp * (ind + 1) + comp] * c;
}

template <typename FType>
void FPlaneAveragingFine<FType>::line_average(
    int comp, amrex::Vector<amrex::Real>& l_vec)
{
    BL_PROFILE("amr-wind::PlaneAveragingFine::line_average");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    for (int i = 0; i < m_ncell_line; i++) {
        l_vec[i] = m_line_average[m_ncomp * i + comp];
    }
}

template <typename FType>
amrex::Real
FPlaneAveragingFine<FType>::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveragingFine::line_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    return m_line_average[m_ncomp * ind + comp];
}

template <typename FType>
void FPlaneAveragingFine<FType>::operator()()
{
    BL_PROFILE("amr-wind::FPlaneAveragingFine::operator");

    m_last_updated_index = m_time.time_index();

    std::fill(m_line_average.begin(), m_line_average.end(), 0.0);

    switch (m_axis) {
    case 0:
        compute_averages(XDir());
        break;
    case 1:
        compute_averages(YDir());
        break;
    case 2:
        compute_averages(ZDir());
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename FType>
template <typename IndexSelector>
void FPlaneAveragingFine<FType>::compute_averages(const IndexSelector& idxOp)
{
    BL_PROFILE("amr-wind::PlaneAveragingFine::compute_averages");

    amrex::AsyncArray<amrex::Real> lavg(
        m_line_average.data(), m_line_average.size());

    amrex::Real* line_avg = lavg.data();
    const int num_comps = m_ncomp;
    const int num_cells = m_ncell_line;
    const amrex::Real line_dx = m_dx;
    const amrex::Real xlo = m_xlo;

    auto g0 = m_field.repo().mesh().Geom(0);
    const amrex::Real denom =
        (amrex::Real)m_ncell_line /
        ((g0.ProbHi(0) - g0.ProbLo(0)) * (g0.ProbHi(1) - g0.ProbLo(1)) *
         (g0.ProbHi(2) - g0.ProbLo(2)));

    const auto& mesh = m_field.repo().mesh();
    const int finestLevel = mesh.finestLevel();

    for (int lev = 0; lev <= finestLevel; ++lev) {

        const auto& geom = m_field.repo().mesh().Geom(lev);
        const amrex::Real dx = geom.CellSize()[m_axis];
        const amrex::Real dy = geom.CellSize()[idxOp.odir1];
        const amrex::Real dz = geom.CellSize()[idxOp.odir2];

        amrex::iMultiFab level_mask;
        if (lev < finestLevel) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& mfab = m_field(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            amrex::Box bx = mfi.tilebox();

            auto fab_arr = mfab.const_array(mfi);
            auto mask_arr = level_mask.const_array(mfi);

            amrex::Box pbx =
                PerpendicularBox<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

            amrex::ParallelFor(
                amrex::Gpu::KernelInfo().setReduction(true), pbx,
                [=] AMREX_GPU_DEVICE(
                    int p_i, int p_j, int p_k,
                    amrex::Gpu::Handler const& handler) noexcept {
                    // Loop over the direction perpendicular to the plane.
                    // This reduces the atomic pressure on the destination
                    // arrays.

                    amrex::Box lbx = ParallelBox<IndexSelector>(
                        bx, amrex::IntVect{p_i, p_j, p_k});

                    for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                        for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                            for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0);
                                 ++i) {

                                if (mask_arr(i, j, k)) {
                                    // cell coordinates
                                    const amrex::Real cell_xlo =
                                        xlo + idxOp(i, j, k) * dx;
                                    const amrex::Real cell_xhi = cell_xlo + dx;

                                    // line indices
                                    const int line_ind_lo = amrex::min(
                                        amrex::max(
                                            static_cast<int>(
                                                (cell_xlo - xlo) / line_dx),
                                            0),
                                        num_cells - 1);
                                    const int line_ind_hi = amrex::min(
                                        amrex::max(
                                            static_cast<int>(
                                                (cell_xhi - xlo) / line_dx),
                                            0),
                                        num_cells - 1);

                                    AMREX_ASSERT(line_ind_lo >= 0);
                                    AMREX_ASSERT(line_ind_hi >= 0);
                                    AMREX_ASSERT(line_ind_lo < num_cells);
                                    AMREX_ASSERT(line_ind_hi < num_cells);

                                    for (int ind = line_ind_lo;
                                         ind <= line_ind_hi; ++ind) {

                                        // line coordinates
                                        const amrex::Real line_xlo =
                                            xlo + ind * line_dx;
                                        const amrex::Real line_xhi =
                                            line_xlo + line_dx;

                                        amrex::Real deltax;

                                        if (line_xlo <= cell_xlo) {
                                            deltax = line_xhi - cell_xlo;
                                        } else if (line_xhi >= cell_xhi) {
                                            deltax = cell_xhi - line_xlo;
                                        } else {
                                            deltax = line_dx;
                                        }
                                        deltax = amrex::min(deltax, dx);
                                        const amrex::Real vol =
                                            deltax * dy * dz;

                                        for (int n = 0; n < num_comps; ++n) {
                                            amrex::Gpu::deviceReduceSum(
                                                &line_avg[num_comps * ind + n],
                                                fab_arr(i, j, k, n) * vol *
                                                    denom,
                                                handler);
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
        }
    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_average.data(), m_line_average.size());
}

template class FPlaneAveragingFine<Field>;
template class FPlaneAveragingFine<ScratchField>;

VelPlaneAveragingFine::VelPlaneAveragingFine(CFDSim& sim, int axis_in)
    : FieldPlaneAveragingFine(
          sim.repo().get_field("velocity"), sim.time(), axis_in),
      temperatureField(sim.repo().get_field("temperature"))
{
    m_line_hvelmag_average.resize(m_ncell_line, 0.0);
    m_line_Su_average.resize(m_ncell_line, 0.0);
    m_line_Sv_average.resize(m_ncell_line, 0.0);
    m_line_Stheta_average.resize(m_ncell_line, 0.0);
}

void VelPlaneAveragingFine::operator()()
{

    BL_PROFILE("amr-wind::VelPlaneAveraging::operator");

    // velocity averages
    FieldPlaneAveragingFine::operator()();

    std::fill(m_line_hvelmag_average.begin(), m_line_hvelmag_average.end(), 0.0);
    std::fill(m_line_Su_average.begin(), m_line_Su_average.end(), 0.0);
    std::fill(m_line_Sv_average.begin(), m_line_Sv_average.end(), 0.0);
    std::fill(m_line_Stheta_average.begin(), m_line_Stheta_average.end(), 0.0);

    switch (m_axis) {
    case 0:
        compute_hvelmag_averages(XDir());
        break;
    case 1:
        compute_hvelmag_averages(YDir());
        break;
    case 2:
        compute_hvelmag_averages(ZDir());
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename IndexSelector>
void VelPlaneAveragingFine::compute_hvelmag_averages(const IndexSelector& idxOp)
{

    BL_PROFILE("amr-wind::VelPlaneAveragingFine::compute_hvelmag_averages");

    amrex::AsyncArray<amrex::Real> lavg_vm(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::AsyncArray<amrex::Real> lavg_Su(
        m_line_Su_average.data(), m_line_Su_average.size());
    amrex::AsyncArray<amrex::Real> lavg_Sv(
        m_line_Sv_average.data(), m_line_Sv_average.size());
    amrex::AsyncArray<amrex::Real> lavg_Stheta(
        m_line_Stheta_average.data(), m_line_Stheta_average.size());

    amrex::Real* line_avg_vm = lavg_vm.data();
    amrex::Real* line_avg_Su = lavg_Su.data();
    amrex::Real* line_avg_Sv = lavg_Sv.data();
    amrex::Real* line_avg_Stheta = lavg_Stheta.data();

    const int num_cells = m_ncell_line;
    const amrex::Real line_dx = m_dx;
    const amrex::Real xlo = m_xlo;

    auto g0 = m_field.repo().mesh().Geom(0);
    const amrex::Real denom =
        (amrex::Real)m_ncell_line /
        ((g0.ProbHi(0) - g0.ProbLo(0)) * (g0.ProbHi(1) - g0.ProbLo(1)) *
         (g0.ProbHi(2) - g0.ProbLo(2)));

    const auto& mesh = m_field.repo().mesh();
    const int finestLevel = mesh.finestLevel();

    for (int lev = 0; lev <= finestLevel; ++lev) {

        const auto& geom = m_field.repo().mesh().Geom(lev);
        const amrex::Real dx = geom.CellSize()[m_axis];
        const amrex::Real dy = geom.CellSize()[idxOp.odir1];
        const amrex::Real dz = geom.CellSize()[idxOp.odir2];

        amrex::iMultiFab level_mask;
        if (lev < finestLevel) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& mfab = m_field(lev);
        const auto& tfab = temperatureField(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            amrex::Box bx = mfi.tilebox();

            auto fab_arr = mfab.const_array(mfi);
            auto t_arr = tfab.const_array(mfi);
            auto mask_arr = level_mask.const_array(mfi);

            amrex::Box pbx =
                PerpendicularBox<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

            amrex::ParallelFor(
                amrex::Gpu::KernelInfo().setReduction(true), pbx,
                [=] AMREX_GPU_DEVICE(
                    int p_i, int p_j, int p_k,
                    amrex::Gpu::Handler const& handler) noexcept {
                    // Loop over the direction perpendicular to the plane.
                    // This reduces the atomic pressure on the destination
                    // arrays.

                    amrex::Box lbx = ParallelBox<IndexSelector>(
                        bx, amrex::IntVect{p_i, p_j, p_k});

                    for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                        for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                            for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0);
                                 ++i) {

                                if (mask_arr(i, j, k)) {
                                    // cell coordinates
                                    const amrex::Real cell_xlo =
                                        xlo + idxOp(i, j, k) * dx;
                                    const amrex::Real cell_xhi = cell_xlo + dx;

                                    // line indices
                                    const int line_ind_lo = amrex::min(
                                        amrex::max(
                                            static_cast<int>(
                                                (cell_xlo - xlo) / line_dx),
                                            0),
                                        num_cells - 1);
                                    const int line_ind_hi = amrex::min(
                                        amrex::max(
                                            static_cast<int>(
                                                (cell_xhi - xlo) / line_dx),
                                            0),
                                        num_cells - 1);

                                    AMREX_ASSERT(line_ind_lo >= 0);
                                    AMREX_ASSERT(line_ind_hi >= 0);
                                    AMREX_ASSERT(line_ind_lo < num_cells);
                                    AMREX_ASSERT(line_ind_hi < num_cells);

                                    for (int ind = line_ind_lo;
                                         ind <= line_ind_hi; ++ind) {

                                        // line coordinates
                                        const amrex::Real line_xlo =
                                            xlo + ind * line_dx;
                                        const amrex::Real line_xhi =
                                            line_xlo + line_dx;

                                        amrex::Real deltax;

                                        if (line_xlo <= cell_xlo) {
                                            deltax = line_xhi - cell_xlo;
                                        } else if (line_xhi >= cell_xhi) {
                                            deltax = cell_xhi - line_xlo;
                                        } else {
                                            deltax = line_dx;
                                        }

                                        deltax = amrex::min(deltax, dx);
                                        const amrex::Real vol =
                                            deltax * dy * dz;

                                        const amrex::Real hvelmag = std::sqrt(
                                            fab_arr(i, j, k, idxOp.odir1) *
                                                fab_arr(i, j, k, idxOp.odir1) +
                                            fab_arr(i, j, k, idxOp.odir2) *
                                                fab_arr(i, j, k, idxOp.odir2));
                                        const amrex::Real Su =
                                            hvelmag *
                                            fab_arr(i, j, k, idxOp.odir1);
                                        const amrex::Real Sv =
                                            hvelmag *
                                            fab_arr(i, j, k, idxOp.odir2);
                                        const amrex::Real Stheta =
                                            hvelmag *
                                            t_arr(i, j, k);

                                        amrex::Gpu::deviceReduceSum(
                                            &line_avg_vm[ind],
                                            hvelmag * vol * denom, handler);
                                        amrex::Gpu::deviceReduceSum(
                                            &line_avg_Su[ind], Su * vol * denom,
                                            handler);
                                        amrex::Gpu::deviceReduceSum(
                                            &line_avg_Sv[ind], Sv * vol * denom,
                                            handler);
                                        amrex::Gpu::deviceReduceSum(
                                            &line_avg_Stheta[ind], Stheta * vol * denom,
                                            handler);
                                    }
                                }
                            }
                        }
                    }
                });
        }
    }

    lavg_vm.copyToHost(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    lavg_Su.copyToHost(m_line_Su_average.data(), m_line_Su_average.size());
    lavg_Sv.copyToHost(m_line_Sv_average.data(), m_line_Sv_average.size());
    lavg_Stheta.copyToHost(m_line_Stheta_average.data(), m_line_Stheta_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_Su_average.data(), m_line_Su_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_Sv_average.data(), m_line_Sv_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_Stheta_average.data(), m_line_Stheta_average.size());
}

amrex::Real
VelPlaneAveragingFine::line_hvelmag_average_interpolated(amrex::Real x) const
{
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_hvelmag_average[ind] * (1.0 - c) +
           m_line_hvelmag_average[ind + 1] * c;
}

amrex::Real
VelPlaneAveragingFine::line_Su_average_interpolated(amrex::Real x) const
{
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_Su_average[ind] * (1.0 - c) + m_line_Su_average[ind + 1] * c;
}

amrex::Real
VelPlaneAveragingFine::line_Sv_average_interpolated(amrex::Real x) const
{
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_Sv_average[ind] * (1.0 - c) + m_line_Sv_average[ind + 1] * c;
}

amrex::Real
VelPlaneAveragingFine::line_Stheta_average_interpolated(amrex::Real x) const
{
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_Stheta_average[ind] * (1.0 - c) + m_line_Stheta_average[ind + 1] * c;
}

} // namespace amr_wind
