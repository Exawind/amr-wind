#include "ThirdMomentAveraging.H"

using namespace amr_wind;

void ThirdMomentAveraging::output_line_average_ascii(
    std::string filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::ThirdMomentAveraging::output_line_average_ascii");

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
                for (int p = 0; p < m_plane_average3.ncomp(); ++p) {
                    outfile << ",<" + m_plane_average1.field().base_name() +
                                   std::to_string(m) + "'" +
                                   m_plane_average2.field().base_name() +
                                   std::to_string(n) + "'" +
                                   m_plane_average3.field().base_name() +
                                   std::to_string(p) + "'>";
                }
            }
        }
        outfile << std::endl;

    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    const int ncomp1 = m_plane_average1.ncomp();
    const int ncomp2 = m_plane_average2.ncomp();
    const int ncomp3 = m_plane_average3.ncomp();

    for (int i = 0; i < m_plane_average1.ncell_line(); ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_plane_average1.line_centroids()[i];
        for (int m = 0; m < ncomp1; ++m)
            for (int n = 0; n < ncomp2; ++n)
                for (int p = 0; p < ncomp3; ++p)
                    outfile << ", " << std::scientific
                            << m_third_moments_line
                                   [m_num_moments * i + ncomp1 * ncomp2 * m +
                                    ncomp2 * n + p];

        outfile << std::endl;
    }
}

void ThirdMomentAveraging::output_line_average_ascii(int step, amrex::Real time)
{
    std::string filename = "third_moment_" + m_plane_average1.field().name() +
                           "_" + m_plane_average2.field().name() + "_" +
                           m_plane_average3.field().name() + ".txt";

    output_line_average_ascii(filename, step, time);
}

ThirdMomentAveraging::ThirdMomentAveraging(
    FieldPlaneAveraging& pa1,
    FieldPlaneAveraging& pa2,
    FieldPlaneAveraging& pa3)
    : m_plane_average1(pa1), m_plane_average2(pa2), m_plane_average3(pa3)
{

    AMREX_ALWAYS_ASSERT(m_plane_average1.axis() == m_plane_average2.axis());
    AMREX_ALWAYS_ASSERT(m_plane_average1.axis() == m_plane_average3.axis());
    AMREX_ALWAYS_ASSERT(m_plane_average1.level() == m_plane_average2.level());
    AMREX_ALWAYS_ASSERT(m_plane_average1.level() == m_plane_average3.level());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_plane() == m_plane_average2.ncell_plane());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_plane() == m_plane_average3.ncell_plane());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_line() == m_plane_average2.ncell_line());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.ncell_line() == m_plane_average3.ncell_line());

    // this could be relaxed if we wanted one plane to be at a different time
    // step than another plane
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.last_updated_index() ==
        m_plane_average2.last_updated_index());
    AMREX_ALWAYS_ASSERT(
        m_plane_average1.last_updated_index() ==
        m_plane_average3.last_updated_index());

    m_num_moments = m_plane_average1.ncomp() * m_plane_average2.ncomp() *
                    m_plane_average3.ncomp();

    m_third_moments_line.resize(
        m_plane_average1.ncell_line() * m_num_moments, 0.0);
}

void ThirdMomentAveraging::operator()()
{

    m_last_updated_index = m_plane_average1.last_updated_index();

    auto& field1 = m_plane_average1.field();
    auto& field2 = m_plane_average2.field();
    auto& field3 = m_plane_average3.field();

    std::fill(m_third_moments_line.begin(), m_third_moments_line.end(), 0.0);

    const int level = m_plane_average1.level();

    switch (m_plane_average1.axis()) {
    case 0:
        compute_average(XDir(), field1(level), field2(level), field3(level));
        break;
    case 1:
        compute_average(YDir(), field1(level), field2(level), field3(level));
        break;
    case 2:
        compute_average(ZDir(), field1(level), field2(level), field3(level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename IndexSelector>
void ThirdMomentAveraging::compute_average(
    const IndexSelector& idxOp,
    const amrex::MultiFab& mfab1,
    const amrex::MultiFab& mfab2,
    const amrex::MultiFab& mfab3)
{

    BL_PROFILE("amr-wind::ThirdMomentAveraging::compute_average");

    amrex::AsyncArray<amrex::Real> lfluc(
        m_third_moments_line.data(), m_third_moments_line.size());
    amrex::Real* line_fluc = lfluc.data();

    amrex::AsyncArray<amrex::Real> lavg1(
        m_plane_average1.line_average().data(),
        m_plane_average1.line_average().size());
    amrex::AsyncArray<amrex::Real> lavg2(
        m_plane_average2.line_average().data(),
        m_plane_average2.line_average().size());
    amrex::AsyncArray<amrex::Real> lavg3(
        m_plane_average3.line_average().data(),
        m_plane_average3.line_average().size());

    amrex::Real* line_avg1 = lavg1.data();
    amrex::Real* line_avg2 = lavg2.data();
    amrex::Real* line_avg3 = lavg3.data();

    amrex::Real denom = 1.0 / (amrex::Real)m_plane_average1.ncell_plane();

    const int ncomp1 = m_plane_average1.ncomp();
    const int ncomp2 = m_plane_average2.ncomp();
    const int ncomp3 = m_plane_average3.ncomp();
    const int nmoments = m_num_moments;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab1, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto mfab_arr1 = mfab1.const_array(mfi);
        auto mfab_arr2 = mfab2.const_array(mfi);
        auto mfab_arr3 = mfab3.const_array(mfi);

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
                        for (int p = 0; p < ncomp3; ++p) {

                            const amrex::Real up3 = mfab_arr3(i, j, k, n) -
                                                    line_avg3[ncomp3 * ind + p];

                            amrex::HostDevice::Atomic::Add(
                                &line_fluc[nmoments * ind + nf],
                                up1 * up2 * up3 * denom);
                            ++nf;
                        }
                    }
                }
            });
    }

    lfluc.copyToHost(m_third_moments_line.data(), m_third_moments_line.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_third_moments_line.data(), m_third_moments_line.size());
}

amrex::Real ThirdMomentAveraging::line_average_interpolated(
    amrex::Real x, int comp1, int comp2, int comp3) const
{
    BL_PROFILE("amr-wind::ThirdMomentAveraging::line_average_interpolated 1");

    AMREX_ALWAYS_ASSERT(comp1 >= 0 && comp1 < m_plane_average1.ncomp());
    AMREX_ALWAYS_ASSERT(comp2 >= 0 && comp2 < m_plane_average2.ncomp());
    AMREX_ALWAYS_ASSERT(comp3 >= 0 && comp3 < m_plane_average3.ncomp());

    const int comp =
        m_plane_average1.ncomp() * m_plane_average2.ncomp() * comp1 +
        m_plane_average2.ncomp() * comp2 + comp3;
    return line_average_interpolated(x, comp);
}

amrex::Real
ThirdMomentAveraging::line_average_interpolated(amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::ThirdMomentAveraging::line_average_interpolated 2");

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

    return m_third_moments_line[m_num_moments * ind + comp] * (1.0 - c) +
           m_third_moments_line[m_num_moments * (ind + 1) + comp] * c;
}

amrex::Real ThirdMomentAveraging::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::ThirdMomentAveraging::line_average_cell 2");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_num_moments);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind + 1 < m_plane_average1.ncell_line());

    return m_third_moments_line[m_num_moments * ind + comp];
}

amrex::Real ThirdMomentAveraging::line_average_cell(
    int ind, int comp1, int comp2, int comp3) const
{
    BL_PROFILE("amr-wind::ThirdMomentAveraging::line_average_cell 1");

    AMREX_ALWAYS_ASSERT(comp1 >= 0 && comp1 < m_plane_average1.ncomp());
    AMREX_ALWAYS_ASSERT(comp2 >= 0 && comp2 < m_plane_average2.ncomp());
    AMREX_ALWAYS_ASSERT(comp3 >= 0 && comp3 < m_plane_average3.ncomp());

    const int comp =
        m_plane_average1.ncomp() * m_plane_average2.ncomp() * comp1 +
        m_plane_average2.ncomp() * comp2 + comp3;

    return line_average_cell(ind, comp);
}
