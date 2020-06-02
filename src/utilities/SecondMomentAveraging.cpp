#include "SecondMomentAveraging.H"

using namespace amr_wind;

void SecondMomentAveraging::output_line_average_ascii(
    std::string filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::SecondMomentAveraging::output_line_average_ascii");

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

        outfile << m_plane_average1.ncell_line() << ", " << m_num_moments + 3
                << std::endl;
        outfile << "#step,time,z";

        for (int m = 0; m < m_plane_average1.ncomp(); ++m) {
            for (int n = 0; n < m_plane_average1.ncomp(); ++n) {
                outfile << ",<" + m_plane_average1.field().base_name() +
                               std::to_string(m) + "'" +
                               m_plane_average2.field().base_name() +
                               std::to_string(n) + "'>";
            }
        }

        outfile << std::endl;

    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    const int ncomp1 = m_plane_average1.ncomp();
    const int ncomp2 = m_plane_average2.ncomp();

    for (int i = 0; i < m_plane_average1.ncell_line(); ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_plane_average1.line_centroids()[i];
        for (int m = 0; m < ncomp1; ++m)
            for (int n = 0; n < ncomp2; ++n)
                outfile << ", " << std::scientific
                        << m_second_moments_line
                               [m_num_moments * i + ncomp1 * m + n];

        outfile << std::endl;
    }
}

void SecondMomentAveraging::output_line_average_ascii(
    int step, amrex::Real time)
{
    std::string filename = "second_moment_" + m_plane_average1.field().name() +
                           "_" + m_plane_average2.field().name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

SecondMomentAveraging::SecondMomentAveraging(
    FieldPlaneAveraging& pa1, FieldPlaneAveraging& pa2)
    : m_plane_average1(pa1), m_plane_average2(pa2)
{

    AMREX_ALWAYS_ASSERT(m_plane_average1.axis() == m_plane_average2.axis());
    AMREX_ALWAYS_ASSERT(m_plane_average1.level() == m_plane_average2.level());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_plane() == m_plane_average2.ncell_plane());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_line() == m_plane_average2.ncell_line());

    // this could be relaxed if we wanted one plane to be at a different time
    // step than another plane
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.last_updated_index() ==
        m_plane_average2.last_updated_index());

    m_num_moments = m_plane_average1.ncomp() * m_plane_average2.ncomp();

    m_second_moments_line.resize(
        m_plane_average1.ncell_line() * m_num_moments, 0.0);
}

void SecondMomentAveraging::operator()()
{

    m_last_updated_index = m_plane_average1.last_updated_index();

    std::fill(m_second_moments_line.begin(), m_second_moments_line.end(), 0.0);

    auto& field1 = m_plane_average1.field();
    auto& field2 = m_plane_average2.field();

    const int level = m_plane_average1.level();

    switch (m_plane_average1.axis()) {
    case 0:
        compute_average(XDir(), field1(level), field2(level));
        break;
    case 1:
        compute_average(YDir(), field1(level), field2(level));
        break;
    case 2:
        compute_average(ZDir(), field1(level), field2(level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename IndexSelector>
void SecondMomentAveraging::compute_average(
    const IndexSelector& idxOp,
    const amrex::MultiFab& mfab1,
    const amrex::MultiFab& mfab2)
{

    BL_PROFILE("amr-wind::SecondMomentAveraging::compute_average");

    amrex::AsyncArray<amrex::Real> lfluc(
        m_second_moments_line.data(), m_second_moments_line.size());
    amrex::Real* line_fluc = lfluc.data();

    amrex::AsyncArray<amrex::Real> lavg1(
        m_plane_average1.line_average().data(),
        m_plane_average1.line_average().size());
    amrex::AsyncArray<amrex::Real> lavg2(
        m_plane_average2.line_average().data(),
        m_plane_average2.line_average().size());

    amrex::Real* line_avg1 = lavg1.data();
    amrex::Real* line_avg2 = lavg2.data();

    amrex::Real denom = 1.0 / (amrex::Real)m_plane_average1.ncell_plane();

    const int ncomp1 = m_plane_average1.ncomp();
    const int ncomp2 = m_plane_average2.ncomp();
    const int nmoments = m_num_moments;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab1, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto mfab_arr1 = mfab1.const_array(mfi);
        auto mfab_arr2 = mfab2.const_array(mfi);

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int ind = idxOp(i, j, k);

                int nf = 0;
                for (int m = 0; m < ncomp1; ++m) {
                    const amrex::Real up1 =
                        mfab_arr1(i, j, k, m) - line_avg1[ncomp1 * ind + m];
                    for (int n = 0; n < ncomp2; ++n) {
                        const amrex::Real up2 =
                            mfab_arr2(i, j, k, n) - line_avg2[ncomp2 * ind + n];
                        amrex::HostDevice::Atomic::Add(
                            &line_fluc[nmoments * ind + nf],
                            up1 * up2 * denom);
                        ++nf;
                    }
                }
            });
    }

    lfluc.copyToHost(
        m_second_moments_line.data(), m_second_moments_line.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_second_moments_line.data(), m_second_moments_line.size());
}

amrex::Real SecondMomentAveraging::line_average_interpolated(
    amrex::Real x, int comp1, int comp2) const
{
    BL_PROFILE("amr-wind::SecondMomentAveraging::line_average_interpolated 1");

    AMREX_ALWAYS_ASSERT(comp1 >= 0 && comp1 < m_plane_average1.ncomp());
    AMREX_ALWAYS_ASSERT(comp2 >= 0 && comp2 < m_plane_average2.ncomp());

    const int comp = m_plane_average1.ncomp() * comp1 + comp2;
    return line_average_interpolated(x, comp);
}

amrex::Real
SecondMomentAveraging::line_average_interpolated(amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::SecondMomentAveraging::line_average_interpolated 2");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_num_moments);

    const amrex::Real dx = m_plane_average1.dx();
    const amrex::Real xlo = m_plane_average1.xlo();
    const int ncell_line = m_plane_average1.ncell_line();

    amrex::Real c = 0.0;
    int ind = 0;

    if (x > xlo + 0.5 * dx) {
        ind = floor((x - xlo) / dx - 0.5);
        const amrex::Real x1 = xlo + (ind + 0.5) * dx;
        c = (x - x1) / dx;
    }

    if (ind + 1 >= ncell_line) {
        ind = ncell_line - 2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 and ind + 1 < ncell_line);

    return m_second_moments_line[m_num_moments * ind + comp] * (1.0 - c) +
           m_second_moments_line[m_num_moments * (ind + 1) + comp] * c;
}

amrex::Real SecondMomentAveraging::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::SecondMomentAveraging::line_average_cell 2");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_num_moments);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind + 1 < m_plane_average1.ncell_line());

    return m_second_moments_line[m_num_moments * ind + comp];
}

amrex::Real
SecondMomentAveraging::line_average_cell(int ind, int comp1, int comp2) const
{
    BL_PROFILE("amr-wind::SecondMomentAveraging::line_average_cell 1");

    AMREX_ALWAYS_ASSERT(comp1 >= 0 && comp1 < m_plane_average1.ncomp());
    AMREX_ALWAYS_ASSERT(comp2 >= 0 && comp2 < m_plane_average2.ncomp());

    const int comp = m_plane_average1.ncomp() * comp1 + comp2;
    return line_average_cell(ind, comp);
}
