//
//  SecondMoment.cpp
//  amr-wind
//

#include "SecondMoment.H"

SecondMoment::SecondMoment(FieldPlaneAveraging& pa1_in, FieldPlaneAveraging& pa2_in)
    : pa1(pa1_in)
    , pa2(pa2_in)
    , axis(pa1_in.get_axis())
    , level(pa1_in.get_level())
{

    //fix me put more checks to make sure these field plane averages are similar, use operator==?
    AMREX_ALWAYS_ASSERT(axis ==  pa2.get_axis());
    AMREX_ALWAYS_ASSERT(level == pa2.get_level());
    AMREX_ALWAYS_ASSERT(pa1.get_ncell_plane() == pa2.get_ncell_plane());
    AMREX_ALWAYS_ASSERT(pa1.get_ncell_line() == pa2.get_ncell_line());

    auto& field1 = pa1.get_field();
    auto& field2 = pa2.get_field();
    ncomp1 = field1.num_comp();
    ncomp2 = field2.num_comp();

    nfluc = ncomp1*ncomp2;

    ncell_line = pa1.get_ncell_line();
    ncell_plane = pa1.get_ncell_plane();

    second_moments_line.resize(ncell_line*nfluc);

    std::fill(second_moments_line.begin(), second_moments_line.end(), 0.0);

    switch (axis) {
        case 0:
            fill_line(XDir(), field1(level), field2(level));
            break;
        case 1:
            fill_line(YDir(), field1(level), field2(level));
            break;
        case 2:
            fill_line(ZDir(), field1(level), field2(level));
            break;
        default:
            amrex::Abort("axis must be equal to 0, 1, or 2");
            break;
    }

}

template<typename IndexSelector>
void SecondMoment::fill_line(const IndexSelector &idxOp, const amrex::MultiFab &mfab1, const amrex::MultiFab &mfab2)
{

    amrex::AsyncArray<amrex::Real> lfluc(second_moments_line.data(), second_moments_line.size());
    amrex::Real *line_fluc = lfluc.data();

    amrex::AsyncArray<amrex::Real> lavg1(pa1.get_line_avg_data(), pa1.get_line_avg_size());
    amrex::AsyncArray<amrex::Real> lavg2(pa2.get_line_avg_data(), pa2.get_line_avg_size());

    amrex::Real *line_avg1 = lavg1.data();
    amrex::Real *line_avg2 = lavg2.data();

    amrex::Real denom = 1.0/(amrex::Real) ncell_plane;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab1, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.tilebox();

        auto mfab_arr1 = mfab1.const_array(mfi);
        auto mfab_arr2 = mfab2.const_array(mfi);

        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ind = idxOp(i,j,k);

            int nf = 0;
            for(int m = 0; m < ncomp1; ++m){
                const amrex::Real up1 = mfab_arr1(i, j, k, m) - line_avg1[ncomp1 * ind + m];
                for(int n = 0; n < ncomp2; ++n){
                    const amrex::Real up2 = mfab_arr2(i, j, k, n) - line_avg2[ncomp2 * ind + n];
                    amrex::HostDevice::Atomic::Add(&line_fluc[nfluc * ind + nf], up1 * up2 *denom);
                    ++nf;
                }
            }

        });

    }

    lfluc.copyToHost(second_moments_line.data(), second_moments_line.size());
    amrex::ParallelDescriptor::ReduceRealSum(second_moments_line.data(), second_moments_line.size());

}

amrex::Real SecondMoment::eval_second_moment(amrex::Real x, int comp1, int comp2) const
{
    BL_PROFILE("amr-wind::SecondMoment::eval_second_moment 1")

    AMREX_ALWAYS_ASSERT(comp1 >= 0 && comp1 < ncomp1);
    AMREX_ALWAYS_ASSERT(comp2 >= 0 && comp2 < ncomp2);

    const int comp = ncomp1*comp1 + comp2;
    return eval_second_moment(x, comp);
}

amrex::Real SecondMoment::eval_second_moment(amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::SecondMoment::eval_second_moment 2")

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < nfluc);

    const amrex::Real dx = pa1.get_dx();
    const amrex::Real xlo = pa1.get_xlo();

    amrex::Real c = 0.0;
    int ind = 0;

    if(x > xlo + 0.5*dx){
        ind = floor((x - xlo)/dx - 0.5);
        const amrex::Real x1 = xlo + (ind + 0.5)*dx;
        c = (x - x1)/dx;
    }

    if( ind+1 >= ncell_line){
        ind = ncell_line - 2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind>=0 and ind+1<ncell_line);

    return second_moments_line[nfluc*ind+comp]*(1.0-c) + second_moments_line[nfluc*(ind+1)+comp]*c;

}
