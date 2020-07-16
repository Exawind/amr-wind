#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include <algorithm>
#include "amr-wind/incflo.H"

namespace amr_wind {

template <typename FType>
void FPlaneAveraging<FType>::output_line_average_ascii(
    std::string filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::FPlaneAveraging::output_line_average_ascii");

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
            outfile << ",<" + m_field.name() + std::to_string(i) + ">";
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
void FPlaneAveraging<FType>::output_line_average_ascii(int step, amrex::Real time)
{
    std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}
    
template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_average_interpolated(amrex::Real x, int comp) const
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
    
template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind < m_ncell_line);

    return m_line_average[m_ncomp * ind + comp];
}

template <typename FType>
void FPlaneAveraging<FType>::compute_line_derivatives()
{
    BL_PROFILE("amr-wind::PlaneAveraging::compute_line_derivatives");
    for (int i = 0; i < m_ncell_line; ++i) {
        for (int n = 0; n < m_ncomp; ++n)
            m_line_deriv[m_ncomp * i + n] =
                line_derivative_of_average_cell(i, n);
    }
}
    
template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_derivative_of_average_cell(
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

template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_derivative_interpolated(
    amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_derivative_interpolated");

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

    return m_line_deriv[m_ncomp * ind + comp] * (1.0 - c) +
        m_line_deriv[m_ncomp * (ind + 1) + comp] * c;
}

template <typename FType>
void FPlaneAveraging<FType>::operator()()
{
    BL_PROFILE("amr-wind::FPlaneAveraging::operator");

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

    if (m_comp_deriv)
        compute_line_derivatives();
    
}

template <typename FType>
template <typename IndexSelector>
void FPlaneAveraging<FType>::compute_averages(
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

void VelPlaneAveraging::operator()()
{
    
    BL_PROFILE("amr-wind::VelPlaneAveraging::operator");

    FieldPlaneAveraging::operator()();

    std::fill(m_line_hvelmag_average.begin(), m_line_hvelmag_average.end(), 0.0);

    switch (m_axis) {
    case 0:
        compute_hvelmag_averages(XDir(), 1, 2, m_field(m_level));
        break;
    case 1:
        compute_hvelmag_averages(YDir(), 0, 2, m_field(m_level));
        break;
    case 2:
        compute_hvelmag_averages(ZDir(), 0, 1, m_field(m_level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }

    if (m_comp_deriv)
        compute_line_hvelmag_derivatives();
    
}

template <typename IndexSelector>
void VelPlaneAveraging::compute_hvelmag_averages(
    const IndexSelector& idx_op,
    const int h1_idx,
    const int h2_idx,
    const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::compute_hvelmag_averages");

    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane;

    amrex::AsyncArray<amrex::Real> lavg(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::Real* line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::ParallelFor(
            bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int ind = idx_op(i, j, k);
                const double hvelmag = amrex::Math::sqrt(
                    fab_arr(i,j,k,h1_idx)*fab_arr(i,j,k,h1_idx) +
                    fab_arr(i,j,k,h2_idx)*fab_arr(i,j,k,h2_idx) );
                amrex::HostDevice::Atomic::Add(
                    &line_avg[ind], hvelmag*denom);
            });
    }

    lavg.copyToHost(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());

}

void VelPlaneAveraging::compute_line_hvelmag_derivatives()
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::compute_line_hvelmag_derivatives");
    for (int i = 0; i < m_ncell_line; ++i) {
            m_line_hvelmag_deriv[i] =
                line_hvelmag_derivative_of_average_cell(i);
    }
}

amrex::Real VelPlaneAveraging::line_hvelmag_derivative_of_average_cell(
    int ind) const
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::line_derivative_of_average_cell");

    AMREX_ALWAYS_ASSERT(ind >= 0 and ind < m_ncell_line);

    amrex::Real dudx;

    if (ind == 0)
        dudx = (m_line_hvelmag_average[(ind + 1)] -
                m_line_hvelmag_average[ind]) /
            m_dx;
    else if (ind == m_ncell_line - 1)
        dudx = (m_line_hvelmag_average[ind] -
                m_line_hvelmag_average[ind - 1]) /
            m_dx;
    else
        dudx = 0.5 *
            (m_line_hvelmag_average[ind + 1] -
             m_line_hvelmag_average[ind - 1]) /
            m_dx;

    return dudx;
}

amrex::Real VelPlaneAveraging::line_hvelmag_average_interpolated(amrex::Real x) const
{

    BL_PROFILE("amr-wind::PlaneAveraging::line_average_interpolated");
    
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

    return m_line_hvelmag_average[ind] * (1.0 - c) +
           m_line_hvelmag_average[ind + 1] * c;
}
    
amrex::Real VelPlaneAveraging::line_hvelmag_average_cell(int ind) const
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::line_hvelmag_average_cell");

    AMREX_ALWAYS_ASSERT(ind >= 0 and ind < m_ncell_line);

    return m_line_average[ind];
}

void VelPlaneAveraging::output_line_average_ascii(
    std::string filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::output_line_average_ascii");

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
        outfile << m_ncell_line << ", " << m_ncomp + 4 << std::endl;
        outfile << "#step,time,z";
        for (int i = 0; i < m_ncomp; ++i)
            outfile << ",<" + m_field.name() + std::to_string(i) + ">";
        outfile << ", <hvelmag>" ;
        outfile << std::endl;
    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    std::cerr << "m_ncomp = " << m_ncomp << std::endl;
    for (int i = 0; i < m_ncell_line; ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_line_xcentroid[i];
        for (int n = 0; n < m_ncomp; ++n) {
            outfile << ", " << std::scientific
                    << m_line_average[m_ncomp * i + n];
        }
        outfile << ", " << std::scientific
                << m_line_hvelmag_average[i];
        outfile << std::endl;
    }
}
    
void VelPlaneAveraging::output_line_average_ascii(int step, amrex::Real time)
{
    std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

}
