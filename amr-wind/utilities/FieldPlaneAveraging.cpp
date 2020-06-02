#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include <algorithm>
#include "amr-wind/incflo.H"

namespace amr_wind {

void FieldPlaneAveraging::output_line_average_ascii(
    std::string filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::FieldPlaneAveraging::output_line_average_ascii");

    if (step != m_last_updated_index) {
        operator()();
    }

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    std::ofstream outfile;
    outfile.precision(m_precision);

    if (step == 1) {
        // make new file
        outfile.open(filename.c_str(), std::ios_base::out);
        outfile << "#ncell,ncomp" << std::endl;
        outfile << m_ncell_line << ", " << m_ncomp + 3 << std::endl;
        outfile << "#step,time,z";
        for (int i = 0; i < m_ncomp; ++i)
            outfile << ",<" + m_field.base_name() + std::to_string(i) + ">";
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

void FieldPlaneAveraging::output_line_average_ascii(int step, amrex::Real time)
{
    std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

amrex::Real FieldPlaneAveraging::line_average_interpolated(amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_average_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    amrex::Real c = 0.0;
    int ind = 0;

    if (x > m_xlo + 0.5 * m_dx) {
        ind = floor((x - m_xlo) / m_dx - 0.5);
        const amrex::Real x1 = m_xlo + (ind + 0.5) * m_dx;
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 and ind + 1 < m_ncell_line);

    return m_line_average[m_ncomp * ind + comp] * (1.0 - c) +
           m_line_average[m_ncomp * (ind + 1) + comp] * c;
}

amrex::Real FieldPlaneAveraging::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind < m_ncell_line);

    return m_line_average[m_ncomp * ind + comp];
}

amrex::Real FieldPlaneAveraging::line_derivative_of_average_cell(
    int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_derivative_of_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind < m_ncell_line);

    amrex::Real dudx;

    if (ind == 0)
        dudx = (m_line_average[m_ncomp * (ind + 1) + comp] -
                m_line_average[m_ncomp * ind + comp]) /
               m_dx;
    else if (ind == m_ncell_line - 1)
        dudx = (m_line_average[m_ncomp * (ind) + comp] -
                m_line_average[m_ncomp * (ind - 1) + comp]) /
               m_dx;
    else
        dudx = 0.5 *
               (m_line_average[m_ncomp * (ind + 1) + comp] -
                m_line_average[m_ncomp * (ind - 1) + comp]) /
               m_dx;

    return dudx;
}

FieldPlaneAveraging::FieldPlaneAveraging(
    amr_wind::Field& field_in, const amr_wind::SimTime& time, int axis_in)
    : m_field(field_in), m_time(time), m_axis(axis_in)
{
    AMREX_ALWAYS_ASSERT(m_axis >= 0 and m_axis <= 2);

    auto geom = m_field.repo().mesh().Geom();

    // level=0 is default, could later make this an input.
    // Might only makes sense for fully covered levels

    m_xlo = geom[m_level].ProbLo(m_axis);
    m_dx = geom[m_level].CellSize(m_axis);

    m_ncomp = m_field.num_comp();

    const amrex::Box& domain = geom[m_level].Domain();
    const amrex::IntVect dom_lo(domain.loVect());
    const amrex::IntVect dom_hi(domain.hiVect());

    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;

    // count number of cells in plane
    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (i != m_axis) m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    m_line_average.resize(m_ncell_line * m_ncomp, 0.0);
    m_line_xcentroid.resize(m_ncell_line);

    for (int i = 0; i < m_ncell_line; ++i) {
        m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dx;
    }
}

void FieldPlaneAveraging::operator()()
{
    BL_PROFILE("amr-wind::FieldPlaneAveraging::operator");

    m_last_updated_index = m_time.time_index();

    std::fill(m_line_average.begin(), m_line_average.end(), 0.0);

    switch (m_axis) {
    case 0:
        compute_averages(XDir(), m_field(m_level));
        break;
    case 1:
        compute_averages(YDir(), m_field(m_level));
        break;
    case 2:
        compute_averages(ZDir(), m_field(m_level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename IndexSelector>
void FieldPlaneAveraging::compute_averages(
    const IndexSelector& idxOp, const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::PlaneAveraging::compute_averages");

    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane;

    amrex::AsyncArray<amrex::Real> lavg(m_line_average.data(), m_line_average.size());

    amrex::Real* line_avg = lavg.data();
    const int ncomp = m_ncomp;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const int ind = idxOp(i, j, k);
                amrex::HostDevice::Atomic::Add(
                    &line_avg[ncomp * ind + n], fab_arr(i, j, k, n)*denom);
            });
    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_average.data(), m_line_average.size());

    // fixme later remove denom above and replace with this
//    std::for_each(m_line_average.begin(), m_line_average.end(), [denom](amrex::Real &el){el *= denom; });

}

}
