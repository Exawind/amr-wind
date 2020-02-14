//
//  PlaneAveraging.cpp
//  incflo
//


#include "PlaneAveraging.H"
#include "incflo.H"

using namespace amrex;

void PlaneAveraging::plot_line_text(std::string filename, int step, Real time)
{
    BL_PROFILE("PlaneAveraging::plot_line_text()");

    if(!ParallelDescriptor::IOProcessor()) return;
    
    std::ofstream outfile;
    outfile.precision(precision);

    if(step == 0){
        // make new file
        outfile.open(filename.c_str(),std::ios_base::out);
        outfile << "# ncell, ncomp" << std::endl;
        outfile << ncell_line << ", " << navg+nfluc+3 << std::endl;
        outfile << "# step, time, z, u_avg, v_avg, w_avg, T_avg, uu, uv, uw, vv, vw, ww, wuu, wuv, wuw, wvv, wvw, www, Tu, Tv, Tw, nu_avg" << std::endl;
    }else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out|std::ios_base::app);
    }

    for(int i=0;i<ncell_line;++i){
        outfile << step << ", " << std::scientific << time << ", " << line_xcentroid[i];
        for(int n=0;n<navg;++n){
            outfile <<  ", " << std::scientific << line_average[navg*i+n];
        }
        
        for(int n=0;n<nfluc;++n){
            outfile <<  ", " << std::scientific << line_fluctuation[nfluc*i+n];
        }
        outfile << std::endl;
    }
}

void PlaneAveraging::plot_line_average_text(std::string filename, int step, Real time)
{
    BL_PROFILE("PlaneAveraging::plot_line_average_text()");

    if(!ParallelDescriptor::IOProcessor()) return;

    std::ofstream outfile;
    outfile.precision(precision);
    
    if(step == 0){
        // make new file
        outfile.open(filename.c_str(),std::ios_base::out);
        outfile << "# ncell, ncomp" << std::endl;
        outfile << ncell_line << ", " << navg+3 << std::endl;
        outfile << "# step, time, z, u_avg, v_avg, w_avg, T_avg, nu_avg" << std::endl;
    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out|std::ios_base::app);
    }

    for(int i=0;i<ncell_line;++i){
        outfile << step << ", " << std::scientific << time << ", " << line_xcentroid[i];
        for(int n=0;n<navg;++n){
            outfile <<  ", " << std::scientific << line_average[navg*i+n];
        }
        outfile << std::endl;
    }
}


void PlaneAveraging::plot_line_binary(std::string filename, int step, Real time)
{
    BL_PROFILE("PlaneAveraging::plot_line_binary()");

    if(!ParallelDescriptor::IOProcessor()) return;

    std::ofstream outfile;
    
    if(step == 0){
        // make new file
        outfile.open(filename.c_str(),std::ios_base::out|std::ios_base::binary);

        outfile.write((char *) &ncell_line, sizeof(int));
        outfile.write((char *) line_xcentroid.data(), sizeof(Real)*line_xcentroid.size());
        const int nc = navg+nfluc+2;
        outfile.write((char *) &nc, sizeof(int));
    } else {
        // append file
        outfile.open(filename.c_str(),std::ios_base::out|std::ios_base::binary|std::ios_base::app);
    }

    outfile.write((char *) &step, sizeof(int));
    outfile.write((char *) &time, sizeof(Real));
    outfile.write((char *) line_average.data(), sizeof(Real)*line_average.size());
    outfile.write((char *) line_fluctuation.data(), sizeof(Real)*line_fluctuation.size());
}

template<typename IndexSelector>
void PlaneAveraging::fill_line(const IndexSelector &idxOp, const amrex::MultiFab& velocity,  const amrex::MultiFab& tracer)
{

    BL_PROFILE("PlaneAveraging::fill_line()");
    
    for(int i=0;i<ncell_line;++i){
        line_xcentroid[i] = xlo + (i+0.5)*dx;
    }
    
    const Real denom = 1.0/(Real) ncell_plane;

    AsyncArray<Real> lavg(line_average.data(), line_average.size());
    AsyncArray<Real> lfluc(line_fluctuation.data(), line_fluctuation.size());
    
    Real *line_average_ = lavg.data();
    Real *line_fluctuation_ = lfluc.data();
    
    int navg_ = navg;
    int nfluc_ = nfluc;
    int u_avg_ = u_avg;
    int v_avg_ = v_avg;
    int w_avg_ = w_avg;
    int T_avg_ = T_avg;
    int uu_ = uu;
    int uv_ = uv;
    int uw_ = uw;
    int vv_ = vv;
    int vw_ = vw;
    int ww_ = ww;
    int wuu_ = wuu;
    int wuv_ = wuv;
    int wuw_ = wuw;
    int wvv_ = wvv;
    int wvw_ = wvw;
    int www_ = www;
    int Tu_ = Tu;
    int Tv_ = Tv;
    int Tw_ = Tw;
 
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(velocity, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        auto vel_arr = velocity.const_array(mfi);
        auto tracer_arr = tracer.const_array(mfi);
//        auto den_arr = density.const_array(mfi);
//        auto eta_arr = eta.const_array(mfi); //fixme eta no longer in global storage, this function could be moved to predictor/corrector functions

        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ind = idxOp(i,j,k);
            HostDevice::Atomic::Add(&line_average_[navg_*ind+u_avg_], vel_arr(i,j,k,0)*denom);
            HostDevice::Atomic::Add(&line_average_[navg_*ind+v_avg_], vel_arr(i,j,k,1)*denom);
            HostDevice::Atomic::Add(&line_average_[navg_*ind+w_avg_], vel_arr(i,j,k,2)*denom);
            HostDevice::Atomic::Add(&line_average_[navg_*ind+T_avg_], tracer_arr(i,j,k,0)*denom);
            // nu+nut = (mu+mut)/rho
//            HostDevice::Atomic::Add(&line_average_[navg_*ind+nu_avg_], eta_arr(i,j,k)/den_arr(i,j,k)*denom);
        });
        
    }
    
    lavg.copyToHost(line_average.data(), line_average.size());
    ParallelDescriptor::ReduceRealSum(line_average.data(), line_average.size());

    AsyncArray<Real> lavg2(line_average.data(), line_average.size());
    line_average_ = lavg2.data();
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(velocity, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        auto vel_arr = velocity.const_array(mfi);
        auto tracer_arr = tracer.const_array(mfi);

        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ind = idxOp(i,j,k);
            // velocity fluctuation
            const Real up = vel_arr(i,j,k,0) - line_average_[navg_*ind+u_avg_];
            const Real vp = vel_arr(i,j,k,1) - line_average_[navg_*ind+v_avg_];
            const Real wp = vel_arr(i,j,k,2) - line_average_[navg_*ind+w_avg_];
            
            //fixme need to enumerate tracer variables too
            // tracer fluctuation
            const Real Tp = tracer_arr(i,j,k,0) - line_average_[navg_*ind+T_avg_];

            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+uu_], up*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+uv_], up*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+uw_], up*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+vv_], vp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+vw_], vp*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+ww_], wp*wp*denom);

            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+wuu_], wp*up*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+wuv_], wp*up*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+wuw_], wp*up*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+wvv_], wp*vp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+wvw_], wp*vp*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+www_], wp*wp*wp*denom);

            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+Tu_], Tp*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+Tv_], Tp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation_[nfluc_*ind+Tw_], Tp*wp*denom);
     
        });
        
    }
    
    lfluc.copyToHost(line_fluctuation.data(), line_fluctuation.size());
    ParallelDescriptor::ReduceRealSum(line_fluctuation.data(), line_fluctuation.size());
        

}

Real PlaneAveraging::eval_line_average(Real x, int comp)
{
    
    Real c = 0.0;
    int ind = 0;
    
    if(x > xlo + 0.5*dx){
        ind = floor((x - xlo)/dx - 0.5);
        const Real x1 = xlo + (ind+0.5)*dx;
        c = (x-x1)/dx;
    }

    AMREX_ALWAYS_ASSERT(ind>=0 and ind<ncell_line);

    return line_average[navg*ind+comp]*(1.0-c) + line_average[navg*(ind+1)+comp]*c;
    
}



PlaneAveraging::PlaneAveraging(amrex::Vector<amrex::Geometry>& geom,
                               amrex::Vector<amrex::MultiFab*> const& velocity,
                               amrex::Vector<amrex::MultiFab*> const& tracer,
                               int axis)
{
    
    AMREX_ALWAYS_ASSERT(axis >=0 and axis <= 2);
    
    // level=0 is default, could later make this an input. Might only makes sense for fully covered levels
    
    xlo = geom[level].ProbLo(axis);
    dx = geom[level].CellSize(axis);
    
    const Box& domain = geom[level].Domain();
    const IntVect dom_lo(domain.loVect());
    const IntVect dom_hi(domain.hiVect());
    
    ncell_line = dom_hi[axis]-dom_lo[axis]+1;

    // count number of cells in plane
    ncell_plane = 1;
    for(int i=0;i<AMREX_SPACEDIM;++i){
        if(i!=axis) ncell_plane *= (dom_hi[i]-dom_lo[i]+1);
    }
    
    line_average.resize(ncell_line*navg,0.0);
    line_fluctuation.resize(ncell_line*nfluc,0.0);
    line_xcentroid.resize(ncell_line,0.0);
    
    switch (axis) {
        case 0:
            fill_line(XDir(),*velocity[level],*tracer[level]);
            break;
        case 1:
            fill_line(YDir(),*velocity[level],*tracer[level]);
            break;
        case 2:
            fill_line(ZDir(),*velocity[level],*tracer[level]);
            break;
        default:
            amrex::Abort("axis must be equal to 0, 1, or 2");
            break;
    }
    
}
