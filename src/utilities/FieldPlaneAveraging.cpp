//
//  FieldPlaneAveraging.cpp
//  amr-wind


#include "FieldPlaneAveraging.H"
#include <algorithm>
#include "incflo.H"

using namespace amrex;
using namespace amr_wind;

Real FieldPlaneAveraging::line_average_interpolated(Real x, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::eval_line_average")

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    Real c = 0.0;
    int ind = 0;

    if(x > m_xlo + 0.5*m_dx){
        ind = floor((x - m_xlo)/m_dx - 0.5);
        const Real x1 = m_xlo + (ind+0.5)*m_dx;
        c = (x-x1)/m_dx;
    }

    if( ind+1 >= m_ncell_line){
        ind = m_ncell_line-2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 and ind + 1 < m_ncell_line);

    return m_line_average[m_ncomp*ind+comp]*(1.0-c) + m_line_average[m_ncomp*(ind+1)+comp]*c;

}

Real FieldPlaneAveraging::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::eval_line_average")

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 and ind+1 < m_ncell_line);

    return m_line_average[m_ncomp*ind+comp];

}


FieldPlaneAveraging::FieldPlaneAveraging(amr_wind::Field &field_in, int axis_in)
    : m_field(field_in)
    , m_axis(axis_in)
{
    AMREX_ALWAYS_ASSERT(m_axis >=0 and m_axis <= 2);

    auto geom = m_field.repo().mesh().Geom();

    // level=0 is default, could later make this an input.
    // Might only makes sense for fully covered levels

    m_xlo = geom[m_level].ProbLo(m_axis);
    m_dx = geom[m_level].CellSize(m_axis);

    m_ncomp = m_field.num_comp();

    const Box& domain = geom[m_level].Domain();
    const IntVect dom_lo(domain.loVect());
    const IntVect dom_hi(domain.hiVect());

    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;

    // count number of cells in plane
    m_ncell_plane = 1;
    for(int i=0;i<AMREX_SPACEDIM;++i){
        if(i!=m_axis) m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    m_line_average.resize(m_ncell_line*m_ncomp);
    m_line_xcentroid.resize(m_ncell_line);

    for(int i=0;i<m_ncell_line;++i){
        m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dx;
    }

}

void FieldPlaneAveraging::operator()()
{
    BL_PROFILE("amr-wind::FieldPlaneAveraging::operator")

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

template<typename IndexSelector>
void FieldPlaneAveraging::compute_averages(const IndexSelector &idxOp,
                                   const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::PlaneAveraging::avg_line")

    const Real denom = 1.0/(Real) m_ncell_plane;

    AsyncArray<Real> lavg(m_line_average.data(), m_line_average.size());

    Real *line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::ParallelFor(bx, m_ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ind = idxOp(i,j,k);
            HostDevice::Atomic::Add(&line_avg[m_ncomp * ind + n], fab_arr(i, j, k, n) * denom);
        });

    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    ParallelDescriptor::ReduceRealSum(m_line_average.data(), m_line_average.size());

}
