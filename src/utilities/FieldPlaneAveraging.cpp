//
//  FieldPlaneAveraging.cpp
//  amr-wind


#include "FieldPlaneAveraging.H"
#include <algorithm>
#include "incflo.H"

using namespace amrex;


Real FieldPlaneAveraging::eval_line_average(Real x, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::eval_line_average")

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < ncomp);

    Real c = 0.0;
    int ind = 0;

    if(x > xlo + 0.5*dx){
        ind = floor((x - xlo)/dx - 0.5);
        const Real x1 = xlo + (ind+0.5)*dx;
        c = (x-x1)/dx;
    }

    if( ind+1 >= ncell_line){
        ind = ncell_line-2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind>=0 and ind+1<ncell_line);

    return line_average[ncomp*ind+comp]*(1.0-c) + line_average[ncomp*(ind+1)+comp]*c;

}

Real FieldPlaneAveraging::eval_line_average(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::eval_line_average")

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < ncomp);
    AMREX_ALWAYS_ASSERT(ind>=0 and ind+1<ncell_line);

    return line_average[ncomp*ind+comp];

}


FieldPlaneAveraging::FieldPlaneAveraging(amr_wind::Field &field_in, int axis_in)
    : field(field_in)
    , axis(axis_in)
{
    AMREX_ALWAYS_ASSERT(axis >=0 and axis <= 2);

    auto geom = field.repo().mesh().Geom();

    // level=0 is default, could later make this an input.
    // Might only makes sense for fully covered levels

    xlo = geom[level].ProbLo(axis);
    dx = geom[level].CellSize(axis);

    ncomp = field.num_comp();

    const Box& domain = geom[level].Domain();
    const IntVect dom_lo(domain.loVect());
    const IntVect dom_hi(domain.hiVect());

    ncell_line = dom_hi[axis] - dom_lo[axis] + 1;

    // count number of cells in plane
    ncell_plane = 1;
    for(int i=0;i<AMREX_SPACEDIM;++i){
        if(i!=axis) ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    line_average.resize(ncell_line*ncomp);
    line_xcentroid.resize(ncell_line);

    for(int i=0;i<ncell_line;++i){
        line_xcentroid[i] = xlo + (i + 0.5) * dx;
    }

}

void FieldPlaneAveraging::operator()()
{
    BL_PROFILE("amr-wind::FieldPlaneAveraging::operator")

    std::fill(line_average.begin(), line_average.end(), 0.0);

    switch (axis) {
        case 0:
            avg_line(XDir(), field(level));
            break;
        case 1:
            avg_line(YDir(), field(level));
            break;
        case 2:
            avg_line(ZDir(), field(level));
            break;
        default:
            amrex::Abort("axis must be equal to 0, 1, or 2");
            break;
    }

}

template<typename IndexSelector>
void FieldPlaneAveraging::avg_line(const IndexSelector &idxOp,
                                   const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::PlaneAveraging::avg_line")

    const Real denom = 1.0/(Real) ncell_plane;

    AsyncArray<Real> lavg(line_average.data(), line_average.size());

    Real *line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int ind = idxOp(i,j,k);
            HostDevice::Atomic::Add(&line_avg[ncomp * ind + n], fab_arr(i, j, k, n) * denom);
        });

    }

    lavg.copyToHost(line_average.data(), line_average.size());
    ParallelDescriptor::ReduceRealSum(line_average.data(), line_average.size());

}
