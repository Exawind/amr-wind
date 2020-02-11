//
//  PlaneAveraging.cpp
//  incflo
//


#include "PlaneAveraging.H"
#include "incflo.H"

using namespace amrex;


int get_index(int axis, int i, int j, int k){

    switch (axis) {
        case 0:
            return i;
        case 1:
            return j;
        case 2:
            return k;
    }
    
}

void PlaneAveraging::plot_line_text(int step, Real time)
{

    std::ofstream outfile;
    std::string filename = "line_plot.txt";
    
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

void PlaneAveraging::plot_line_average_text(int step, Real time)
{

    std::ofstream outfile;
    std::string filename = "line_average_plot.txt";

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


void PlaneAveraging::plot_line_binary(int step, Real time)
{

    std::ofstream outfile;

    std::string filename = "line_plot.bin";
    
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

void PlaneAveraging::fill_line(const amrex::MultiFab& velocity,  const amrex::MultiFab& tracer)
{

    for(int i=0;i<ncell_line;++i){
        line_xcentroid[i] = xlo + (i+0.5)*dx;
    }
    
    const Real denom = 1.0/(Real) ncell_plane;
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(velocity, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        auto vel_arr = velocity.const_array(mfi);
        auto tracer_arr = tracer.const_array(mfi);
//        const auto& eta_arr = ld.eta.array(mfi); //fixme eta no longer in global storage, this function could be moved to predictor/corrector functions

        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ind = get_index(axis,i,j,k);
            HostDevice::Atomic::Add(&line_average[navg*ind+u_avg], vel_arr(i,j,k,0)*denom);
            HostDevice::Atomic::Add(&line_average[navg*ind+v_avg], vel_arr(i,j,k,1)*denom);
            HostDevice::Atomic::Add(&line_average[navg*ind+w_avg], vel_arr(i,j,k,2)*denom);
            HostDevice::Atomic::Add(&line_average[navg*ind+T_avg], tracer_arr(i,j,k,0)*denom);
            
        });
        
    }
    
    ParallelDescriptor::ReduceRealSum(line_average.data(), line_average.size());

    
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
            
            const int ind = get_index(axis,i,j,k);
        
            const Real up = vel_arr(i,j,k,0) - line_average[navg*ind+u_avg];
            const Real vp = vel_arr(i,j,k,1) - line_average[navg*ind+v_avg];
            const Real wp = vel_arr(i,j,k,2) - line_average[navg*ind+w_avg];
            const Real Tp = tracer_arr(i,j,k,0) - line_average[navg*ind+T_avg];

            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+uu], up*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+uv], up*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+uw], up*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+vv], vp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+vw], vp*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+ww], wp*wp*denom);

            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+wuu], wp*up*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+wuv], wp*up*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+wuw], wp*up*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+wvv], wp*vp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+wvw], wp*vp*wp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+www], wp*wp*wp*denom);

            
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+Tu], Tp*up*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+Tv], Tp*vp*denom);
            HostDevice::Atomic::Add(&line_fluctuation[nfluc*ind+Tw], Tp*wp*denom);

            // nu+nut = (mu+mut)/rho
//            mfab_arr(i,j,k,nu_avg) = eta_arr(i,j,k)/den_arr(i,j,k);
     
        });
        
    }
    
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

PlaneAveraging::PlaneAveraging(incflo& a_incflo,
                               amrex::Vector<amrex::MultiFab*> const& velocity,
                               amrex::Vector<amrex::MultiFab*> const& tracer,
                               int axis)
    : m_incflo(a_incflo), axis(axis)
{
    
    AMREX_ALWAYS_ASSERT(axis >=0 and axis <= 2);
    
    xlo = m_incflo.Geom(level).ProbLo(axis);
    dx = m_incflo.Geom(level).CellSize(axis);
    
    const Box& domain = m_incflo.Geom(level).Domain();
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

    fill_line(*velocity[level],*tracer[level]);
    
}


//void PlaneAveraging::single_level_line_plot(int s, int b, int axis, int ncomp, Vector<Real> &line){
//
//    IntVect small(0,0,0);
//    IntVect big(0,0,0);
//
//    small.setVal(axis,s);
//    big.setVal(axis,b);
//
//    // have every processor make a 1d box
//    Box bx1d = Box(small,big);
//
//     // make a box array with a single box
//    BoxArray ba1d(bx1d);
//
//    DistributionMapping dm1d;// {ba1d}; // could use the boxarray but not sure if always guaranteed that proc=0 owns it
//    // only processor 0 owns the box
//    Vector<int> pmap {0};
//    dm1d.define(pmap);
//
//    const int ngrow = 0;
//
//     // create a multiFAB using one box on one proc
//    MultiFab mfab1d(ba1d,dm1d,ncomp,ngrow);
//
//    switch(axis){
//        case 0:
//            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//            {
//                Box bx = mfi.tilebox();
//                auto fab_arr = mfab1d.array(mfi);
//                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
//                {
//                    fab_arr(i,0,0,n) = line[n + ncomp*i];
//                });
//            }
//            break;
//        case 1:
//            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//            {
//                Box bx = mfi.tilebox();
//                auto fab_arr = mfab1d.array(mfi);
//                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
//                {
//                    fab_arr(0,j,0,n) = line[n + ncomp*j];
//                });
//            }
//            break;
//        case 2:
//            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
//            {
//                Box bx = mfi.tilebox();
//                auto fab_arr = mfab1d.array(mfi);
//                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
//                {
//                    fab_arr(0,0,k,n) = line[n + ncomp*k];
//                });
//            }
//            break;
//    }
//
//    std::string line_plot_file{"line_plot"};
//
//    const std::string& plotfilename = amrex::Concatenate(line_plot_file, m_nstep);
//
//    amrex::Print() << "  Writing line plot file " << plotfilename << " at time " << m_cur_time << std::endl;
//
//    Vector<std::string> pltscaVarsName;
//
//    pltscaVarsName.push_back("<u>");
//    pltscaVarsName.push_back("<v>");
//    pltscaVarsName.push_back("<w>");
//    pltscaVarsName.push_back("<T>");
//
//    pltscaVarsName.push_back("<u'u'>");
//    pltscaVarsName.push_back("<u'v'>");
//    pltscaVarsName.push_back("<u'w'>");
//    pltscaVarsName.push_back("<v'v'>");
//    pltscaVarsName.push_back("<v'w'>");
//    pltscaVarsName.push_back("<w'w'>");
//
//    pltscaVarsName.push_back("<w'u'u'>");
//    pltscaVarsName.push_back("<w'u'v'>");
//    pltscaVarsName.push_back("<w'u'w'>");
//    pltscaVarsName.push_back("<w'v'v'>");
//    pltscaVarsName.push_back("<w'v'w'>");
//    pltscaVarsName.push_back("<w'w'w'>");
//
//    pltscaVarsName.push_back("<T'u'>");
//    pltscaVarsName.push_back("<T'v'>");
//    pltscaVarsName.push_back("<T'w'>");
//    pltscaVarsName.push_back("<nu+nu_t>");
//
//    //fixme todo add Rij, qj, nu_SGS https://a2ehfm-ecp.slack.com/archives/C3V26K34G/p1522245519000450
//
//    int level_step = 0;
//
//    // in debug mode amrex output will complain if these are not the same size
//    AMREX_ALWAYS_ASSERT(mfab1d.nComp() == pltscaVarsName.size());
//
//    amrex::WriteSingleLevelPlotfile(plotfilename, mfab1d, pltscaVarsName, geom[0], m_cur_time, level_step);
//
//    WriteJobInfo(plotfilename);
//
//}
