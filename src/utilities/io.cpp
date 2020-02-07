#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include <incflo.H>

using namespace amrex;

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void incflo::WriteHeader(const std::string& name, bool is_checkpoint) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        std::string HeaderFileName(name + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(),
                        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if(!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        if(is_checkpoint) {
            HeaderFile << "Checkpoint version: 1\n";
        } else {
            HeaderFile << "HyperCLaw-V1.1\n";
        }

        HeaderFile << finest_level << "\n";

        // Time stepping controls
        HeaderFile << m_nstep << "\n";
        HeaderFile << m_cur_time << "\n";
        HeaderFile << m_dt << "\n";
        HeaderFile << m_prev_dt << "\n";
        HeaderFile << m_prev_prev_dt << "\n";

        // Geometry
        for(int i = 0; i < BL_SPACEDIM; ++i) {
            HeaderFile << Geom(0).ProbLo(i) << ' ';
        }
        HeaderFile << '\n';

        for(int i = 0; i < BL_SPACEDIM; ++i)
            HeaderFile << Geom(0).ProbHi(i) << ' ';
        HeaderFile << '\n';

        // BoxArray
        for(int lev = 0; lev <= finest_level; ++lev)
        {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }
}

void incflo::WriteCheckPointFile() const
{
    BL_PROFILE("incflo::WriteCheckPointFile()");

    const std::string& checkpointname = amrex::Concatenate(m_check_file, m_nstep);

    amrex::Print() << "\n\t Writing checkpoint " << checkpointname << std::endl;

    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

    bool is_checkpoint = true;
    WriteHeader(checkpointname, is_checkpoint);
    WriteJobInfo(checkpointname);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Write(m_leveldata[lev]->velocity,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "velocity"));

        VisMF::Write(m_leveldata[lev]->density,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "density"));

        if (m_ntrac > 0) {
            VisMF::Write(m_leveldata[lev]->tracer,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "tracer"));
        }

        VisMF::Write(m_leveldata[lev]->gp,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gradp"));

        VisMF::Write(m_leveldata[lev]->p,
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p"));
    }
}

void incflo::ReadCheckpointFile()
{
    BL_PROFILE("incflo::ReadCheckpointFile()");

    amrex::Print() << "Restarting from checkpoint " << m_restart_file << std::endl;

    Real prob_lo[BL_SPACEDIM];
    Real prob_hi[BL_SPACEDIM];

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(m_restart_file + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file 
    
    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    GotoNextLine(is);

    // Step count
    is >> m_nstep;
    GotoNextLine(is);

    // Current time
    is >> m_cur_time;
    GotoNextLine(is);

    // Time step size
    is >> m_dt;
    GotoNextLine(is);

    is >> m_prev_dt;
    GotoNextLine(is);

    is >> m_prev_prev_dt;
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while(lis >> word)
        {
            prob_hi[i++] = std::stod(word);
        }
    }

    // Set up problem domain
    RealBox rb(prob_lo, prob_hi);
    Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(lev, Geometry(Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                                  Geom(lev).isPeriodic()));
    }

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // Create distribution mapping
        DistributionMapping dm{ba, ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, m_cur_time, ba, dm);
    }

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

    // Load the field data
    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Read(m_leveldata[lev]->velocity,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "velocity"));

        VisMF::Read(m_leveldata[lev]->density,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "density"));

        if (m_ntrac > 0) {
            VisMF::Read(m_leveldata[lev]->tracer,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "tracer"));
        }

        VisMF::Read(m_leveldata[lev]->gp,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "gradp"));

        VisMF::Read(m_leveldata[lev]->p,
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "p"));
    }

    amrex::Print() << "Restart complete" << std::endl;
}

void incflo::WriteJobInfo(const std::string& path) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = path;
        std::string PrettyLine =
            "===============================================================================\n";

        FullPathJobInfoFile += "/incflo_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " incflo Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

        jobInfoFile << "\n\n";

        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
        jobInfoFile << "FCOMP:         " << buildInfoGetFcomp() << "\n";
        jobInfoFile << "FCOMP version: " << buildInfoGetFcompVersion() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = buildInfoGetGitHash(1);
        const char* githash2 = buildInfoGetGitHash(2);
        if(std::strlen(githash1) > 0)
        {
            jobInfoFile << "incflo git hash: " << githash1 << "\n";
        }
        if(std::strlen(githash2) > 0)
        {
            jobInfoFile << "AMReX git hash: " << githash2 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for(int i = 0; i <= finest_level; i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for(int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                jobInfoFile << geom[i].Domain().length(dir) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << "\n\n";

        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}

void incflo::WritePlotFile() 
{
    BL_PROFILE("incflo::WritePlotFile()");

    if (m_plt_vort or m_plt_divu) {
        for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_EB
            const int ng = (EBFactory(0).isAllRegular()) ? 1 : 2;
#else
            const int ng = 1;
#endif
            fillpatch_velocity(lev, m_cur_time, m_leveldata[lev]->velocity, ng);
        }
    }


    const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

    amrex::Print() << "  Writing plotfile " << plotfilename << " at time " << m_cur_time << std::endl;

    int ncomp = 0;

    // Velocity components
    if (m_plt_velx) ++ncomp;
    if (m_plt_vely) ++ncomp;
    if (m_plt_velz) ++ncomp;
    // Pressure gradient components
    if (m_plt_gpx) ++ncomp;
    if (m_plt_gpy) ++ncomp;
    if (m_plt_gpz) ++ncomp;
    // Density
    if (m_plt_rho) ++ncomp;
    // Tracers
    if (m_plt_tracer) ncomp += m_ntrac;
    // Pressure
    if(m_plt_p) ++ncomp;
    // Apparent viscosity
    if(m_plt_eta) ++ncomp;
    // Vorticity
    if(m_plt_vort) ++ncomp;
    // Magnitude of the rate-of-strain tensor 
    if(m_plt_strainrate) ++ncomp;
    // Magnitude of the stress tensor 
    if(m_plt_stress) ++ncomp;
    // Divergence of velocity field
    if(m_plt_divu) ++ncomp;
#ifdef AMREX_USE_EB
    // Cut cell volume fraction
    if(m_plt_vfrac) ++ncomp;
#endif

    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
    }

    Vector<std::string> pltscaVarsName;
    int icomp = 0;
    if (m_plt_velx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velx");
        ++icomp;
    }
    if (m_plt_vely) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vely");
        ++icomp;
    }
    if (m_plt_velz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->velocity, 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velz");
        ++icomp;
    }
    if (m_plt_gpx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpx");
        ++icomp;
    }
    if (m_plt_gpy) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpy");
        ++icomp;
    }
    if (m_plt_gpz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->gp, 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpz");
        ++icomp;
    }
    if (m_plt_rho) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->density, 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("density");
        ++icomp;
    }
    if (m_plt_tracer) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], m_leveldata[lev]->tracer, 0, icomp, m_ntrac, 0);
        }
        for (int i = 0; i < m_ntrac; ++i) {
            pltscaVarsName.push_back("tracer"+std::to_string(i));
        }
        icomp += m_ntrac;
    }
    if (m_plt_p) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::average_node_to_cellcenter(mf[lev], icomp, m_leveldata[lev]->p, 0, 1);
        }
        pltscaVarsName.push_back("p");
        ++icomp;
    }
    if (m_plt_eta) {
        amrex::Abort("plt_eta: xxxxx TODO");
        pltscaVarsName.push_back("eta");
        ++icomp;
    }
    if (m_plt_vort) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab vort(mf[lev], amrex::make_alias, icomp, 1);
            ComputeVorticity(lev, m_cur_time, vort, m_leveldata[lev]->velocity);
        }
        pltscaVarsName.push_back("vort");
        ++icomp;
    }
    if (m_plt_strainrate) {
        amrex::Abort("plt_strainrate: xxxxx TODO");
        pltscaVarsName.push_back("strainrate");
        ++icomp;
    }
    if (m_plt_stress) {
        amrex::Abort("plt_stress: xxxxx TODO");
        pltscaVarsName.push_back("stress");
        ++icomp;
    }
    if (m_plt_divu) {
        amrex::Abort("plt_divu: xxxxx TODO");
        pltscaVarsName.push_back("divu");
        ++icomp;
    }
#ifdef AMREX_USE_EB
    if (m_plt_vfrac) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], EBFactory(lev).getVolFrac(), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vfrac");
        ++icomp;
    }
#endif

#ifdef AMREX_USE_EB
    for (int lev = 0; lev <= finest_level; ++lev) {
        EB_set_covered(mf[lev], 0.0);
    }
#endif

    AMREX_ALWAYS_ASSERT(ncomp == static_cast<int>(pltscaVarsName.size()));

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
    // but will never change unless we use subcycling. 
    // If we do use subcycling, this should be a incflo class member. 
    Vector<int> istep(finest_level + 1, m_nstep);

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf), 
                                   pltscaVarsName, Geom(), m_cur_time, istep, refRatio());
    WriteJobInfo(plotfilename);
}


enum avg_var {u_avg, v_avg, w_avg, T_avg, uu, uv, uw, vv, vw, ww, wuu, wuv, wuw, wvv, wvw, www, Tu, Tv, Tw, nu_avg, last_avg_var=nu_avg};

void incflo::set_mfab_spatial_averaging_quantities(MultiFab &mfab, int lev, int axis, Vector<Real> &line)
{

    BL_PROFILE("incflo::set_mfab_spatial_averaging_quantities()");
    
    const int ncomp = mfab.nComp();
    
    auto& ld = *m_leveldata[lev];
               
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();

        const auto& mfab_arr = mfab.array(mfi);
        const auto& vel_arr = ld.velocity.array(mfi);
        const auto& tracer_arr = ld.tracer.array(mfi);
//        const auto& eta_arr = ld.eta.array(mfi); //fixme eta no longer in global storage, this function could be moved to predictor/corrector functions
        const auto& den_arr = ld.density.array(mfi);

         // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
        amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // velocities
            mfab_arr(i,j,k,u_avg) = vel_arr(i,j,k,0);
            mfab_arr(i,j,k,v_avg) = vel_arr(i,j,k,1);
            mfab_arr(i,j,k,w_avg) = vel_arr(i,j,k,2);

            // potential temperature
            mfab_arr(i,j,k,T_avg) = tracer_arr(i,j,k,0);
            
            int ind;
            if(axis==0){
                ind = i;
            } else if(axis==1){
                ind = j;
            } else {
                ind = k;
            }
            
            // fluctuations
            const Real up = vel_arr(i,j,k,0) - line[ncomp*ind+u_avg];
            const Real vp = vel_arr(i,j,k,1) - line[ncomp*ind+v_avg];
            const Real wp = vel_arr(i,j,k,2) - line[ncomp*ind+w_avg];
            const Real Tp = tracer_arr(i,j,k,0) - line[ncomp*ind+T_avg];
            
            mfab_arr(i,j,k,uu) = up*up;
            mfab_arr(i,j,k,uv) = up*vp;
            mfab_arr(i,j,k,uw) = up*wp;
            mfab_arr(i,j,k,vv) = vp*vp;
            mfab_arr(i,j,k,vw) = vp*wp;
            mfab_arr(i,j,k,ww) = wp*wp;
            
            mfab_arr(i,j,k,wuu) = wp*up*up;
            mfab_arr(i,j,k,wuv) = wp*up*vp;
            mfab_arr(i,j,k,wuw) = wp*up*wp;
            mfab_arr(i,j,k,wvv) = wp*vp*vp;
            mfab_arr(i,j,k,wvw) = wp*vp*wp;
            mfab_arr(i,j,k,www) = wp*wp*wp;
            
            mfab_arr(i,j,k,Tu) = Tp*up;
            mfab_arr(i,j,k,Tv) = Tp*vp;
            mfab_arr(i,j,k,Tw) = Tp*wp;
            
            // nu+nut = (mu+mut)/rho
//            mfab_arr(i,j,k,nu_avg) = eta_arr(i,j,k)/den_arr(i,j,k);

         });
    }
    
}


void incflo::plane_average(const MultiFab& mfab, const Real nplane_cells, const int axis, const int ncomp, Vector<Real> &line){

    BL_PROFILE("incflo::plane_average()");
  
    switch (axis) {
        case 0:
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                const auto mfab_arr = mfab.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n, {HostDevice::Atomic::Add(&line[ncomp*i+n], mfab_arr(i,j,k,n)/nplane_cells);} );
            }

            break;

        case 1:
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                const auto mfab_arr = mfab.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n, {HostDevice::Atomic::Add(&line[ncomp*j+n], mfab_arr(i,j,k,n)/nplane_cells);} );
            }

            break;

        case 2:
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                const auto mfab_arr = mfab.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n, {HostDevice::Atomic::Add(&line[ncomp*k+n], mfab_arr(i,j,k,n)/nplane_cells);} );
            }

            break;

        default:
            amrex::Abort("aborting axis should be between 0 and 2 \n");
            break;
    }

}


void incflo::set_average_quantities(const int s, const int b, const int axis, const int ncomp, amrex::Vector<Real> &line)
{

    AMREX_ALWAYS_ASSERT(s>=0);
    AMREX_ALWAYS_ASSERT(b>s);
  
    const Real dz = geom[0].CellSize(axis);
    const Real zlo = geom[0].ProbLo(axis);
    
    const Real half_height = zlo + 0.5*dz;

    vx_mean_ground = line[ncomp*s+u_avg];
    vy_mean_ground = line[ncomp*s+v_avg];
    
//    nu_mean_ground = line[ncomp*s+nu_avg];
//    // fixme circular dependency so need to hack this for now
//    if(nstep == -1) nu_mean_ground = 1.0;
//    amrex::Print() << "nu mean ground: " << nu_mean_ground << std::endl;
    
    amrex::Print() << "ground half cell height, vx, vy: " << half_height << ' ' << vx_mean_ground  << ' ' << vy_mean_ground << std::endl;
   
    // search for cells near log_law_sampling_height
    Real height = half_height;
    Real c = 0.0;
    int ind = s;
    
    if(m_log_law_sampling_height > half_height){
        height = m_log_law_sampling_height;
        ind = s + floor((m_log_law_sampling_height - zlo)/dz - 0.5);
        const Real z1 = zlo + (ind+0.5)*dz;
        c = (m_log_law_sampling_height-z1)/dz;
    }

    if(ind   < 0) amrex::Abort("low_law_sampling_height is set incorrectly since conversion to integer is negative \n");
    if(ind+1 > b) amrex::Abort("low_law_sampling_height is set incorrectly since conversion to integer is larger than domain \n");

    const Real vx = line[ncomp*ind+u_avg]*(1.0-c) + line[ncomp*(ind+1)+u_avg]*c;
    const Real vy = line[ncomp*ind+v_avg]*(1.0-c) + line[ncomp*(ind+1)+v_avg]*c;

    const Real uh = sqrt(pow(vx,2) + pow(vy,2));
        
    // simple shear stress model for neutral BL
    // apply as an inhomogeneous Neumann BC
    utau = m_kappa*uh/log(height/m_surface_roughness_z0);
    
    amrex::Print() << "log law sampling height, u, utau: " << height << ' ' << uh  << ' ' << utau << std::endl;
      
    // search for cells near abl_forcing_height
    c = 0.0;
    ind = s;
    if(m_abl_forcing_height > half_height){
        height = m_abl_forcing_height;
        ind = s + floor((m_abl_forcing_height - zlo)/dz - 0.5);
        const Real z1 = zlo + (ind+0.5)*dz;
        c = (m_abl_forcing_height-z1)/dz;
    }

    if(ind   < 0) amrex::Abort("abl_forcing_height is set incorrectly since conversion to integer is negative \n");
    if(ind+1 > b) amrex::Abort("abl_forcing_height is set incorrectly since conversion to integer is larger than domain \n");
   
    vx_mean = line[ncomp*ind+u_avg]*(1.0-c) + line[ncomp*(ind+1)+u_avg]*c;
    vy_mean = line[ncomp*ind+v_avg]*(1.0-c) + line[ncomp*(ind+1)+v_avg]*c;

    
    amrex::Print() << "abl forcing height z, vx, vy: " << height << ' ' << vx_mean << ' ' << vy_mean << std::endl;
    
}


void incflo::spatially_average_quantities_down(int plot_type)
{

    BL_PROFILE("incflo::spatially_average_quantities_down()");

    // for now only have this set to level 0
    const int lev = 0;
    // fixme this could be an input if useful for other problems
    const int axis = 2;
    
    // must match number of quantities in set_mfab_spatial_averaging_quantities(...)
    const int ncomp = last_avg_var+1;

    const Box& domain = geom[lev].Domain();
    const IntVect dom_lo(domain.loVect());
    const IntVect dom_hi(domain.hiVect());
    
    const int s = domain.smallEnd(axis);
    const int b = domain.bigEnd(axis);
    const int ncell = b-s+1;
    
    // code below depends on small index s starting at 0, it would be weird if it didn't
    AMREX_ALWAYS_ASSERT(s==0);
    
    amrex::Vector<Real> line(ncell*ncomp,0.0);

    // count number of cells in plane
    Real nplane_cells = 1.0;
    for(int i=0;i<AMREX_SPACEDIM;++i){
        if(i!=axis) nplane_cells *= (dom_hi[i]-dom_lo[i]+1);
    }
         
    // no ghost cells needed
    const int ngrow = 0;
    
    // first pass is to get averages
    MultiFab mfab(grids[lev], dmap[lev], ncomp, ngrow);
    set_mfab_spatial_averaging_quantities(mfab,lev,axis,line);
    // ncomp is used here and is safer but it could actually be 4 depending on the ordering of average down function
    plane_average(mfab, nplane_cells, axis, ncomp, line);
    // sum all lines together
    ParallelDescriptor::ReduceRealSum(line.data(), line.size());
   
    // second pass has averages and can now calculate fluctuations
    set_mfab_spatial_averaging_quantities(mfab, lev, axis, line);
    // reset back to zero
    line = Vector<Real> (line.size(),0.0);
    plane_average(mfab, nplane_cells, axis, ncomp, line);
    // sum all lines together
    ParallelDescriptor::ReduceRealSum(line.data(), line.size());

    // set all the stored average quantities for every time step
    // fixme warning this has specific ABL stuff in it
    if(m_probtype==35) set_average_quantities(s, b, axis, ncomp, line);
    
    if(!(m_line_plot_int > 0 and (m_nstep % m_line_plot_int == 0))) return;

    // single text file for line data
    // fixme does this deserve it's own function?
    if(plot_type==0 and ParallelDescriptor::IOProcessor()){
        std::ofstream outfile;
        outfile.precision(4);//fixme how much precision do we want?
        if(m_nstep == 0){
            // make new file
            outfile.open("line_plot.txt",std::ios_base::out);
            outfile << "# ncell, ncomp" << std::endl;
            outfile << ncell << ", " << ncomp << std::endl;
            outfile << "# step, time, z, u_avg, v_avg, w_avg, T_avg, uu, uv, uw, vv, vw, ww, wuu, wuv, wuw, wvv, wvw, www, Tu, Tv, Tw, nu_avg" << std::endl;
        }else {
            // append file
            outfile.open("line_plot.txt", std::ios_base::out|std::ios_base::app);
        }
        
        auto dx = geom[lev].CellSizeArray();

        for(int i=s;i<=b;++i){
            Real z = (i+0.5)*dx[axis];
            outfile << m_nstep << ", " << std::scientific << m_cur_time << ", " << z;
            for(int n=0;n<ncomp;++n){
                outfile <<  ", " << std::scientific << line[n + ncomp*i];
            }
            outfile << std::endl;
        }

    }
    
    
    // single binary file for line data
    if(plot_type == 1 and ParallelDescriptor::IOProcessor()){
        std::ofstream outfile;

        if(m_nstep == 0){
            
            // make new file
            outfile.open("line_plot.bin",std::ios_base::out|std::ios_base::binary);
            
            outfile.write((char *) &ncell, sizeof(int));
            outfile.write((char *) &ncomp, sizeof(int));
            
            auto dx = geom[lev].CellSizeArray();

            for(int i=s;i<=b;++i){
                Real z = (i+0.5)*dx[axis];
                outfile.write((char *) &z, sizeof(Real));
            }
            
        }else {
            outfile.open("line_plot.bin",std::ios_base::out|std::ios_base::binary|std::ios_base::app);
        }

        outfile.write((char *) &m_nstep, sizeof(int));
        outfile.write((char *) &m_cur_time, sizeof(Real));
        outfile.write((char *) line.data(), sizeof(Real)*ncell*ncomp);
    }
    
    
    // binary file test
    if(plot_type == 1 and ParallelDescriptor::IOProcessor()){
        std::ifstream infile;

        infile.open("line_plot.bin",std::ios_base::in|std::ios_base::binary);
        if(!infile) amrex::Abort("line plot binary file not found");

        int ncell_in;
        int ncomp_in;
        infile.read((char *) &ncell_in, sizeof(int));
        infile.read((char *) &ncomp_in, sizeof(int));
        
        AMREX_ALWAYS_ASSERT(ncell_in == b-s+1);
        AMREX_ALWAYS_ASSERT(ncomp_in == ncomp);

        Real z[ncell_in];
        infile.read((char *) z, sizeof(Real)*ncell_in);
        auto dx = geom[lev].CellSizeArray();
        for(int i=s,j=0;i<=b;++i,++j)
        {
            Real zdiff = fabs((i+0.5)*dx[axis]- z[j]);
            if(zdiff > 1.0e-12) amrex::Abort("z diff failed");
        }
        
        int nstep_in;
        Real cur_time_in;
        
        infile.read((char *) &nstep_in, sizeof(int));
        infile.read((char *) &cur_time_in, sizeof(Real));
        
        if(m_nstep==0){
            Real line_in[ncell_in*ncomp_in];
            infile.read((char *) line_in, sizeof(Real)*ncell*ncomp);
            for(int i=0;i<ncell;++i){
                for(int n=0;n<ncomp;++n){
                    Real ldiff = fabs(line[n + ncomp*i]-line_in[n+ncomp*i]);
                    if(ldiff > 1.0e-12) amrex::Abort("line diff failed");

                }
            }
        }        
    }
    
    if(plot_type==2) single_level_line_plot(s, b, axis, ncomp, line);
    

 }

void incflo::single_level_line_plot(int s, int b, int axis, int ncomp, Vector<Real> &line){
    
    IntVect small(0,0,0);
    IntVect big(0,0,0);

    small.setVal(axis,s);
    big.setVal(axis,b);

    // have every processor make a 1d box
    Box bx1d = Box(small,big);
            
     // make a box array with a single box
    BoxArray ba1d(bx1d);

    DistributionMapping dm1d;// {ba1d}; // could use the boxarray but not sure if always guaranteed that proc=0 owns it
    // only processor 0 owns the box
    Vector<int> pmap {0};
    dm1d.define(pmap);

    const int ngrow = 0;
    
     // create a multiFAB using one box on one proc
    MultiFab mfab1d(ba1d,dm1d,ncomp,ngrow);

    switch(axis){
        case 0:
            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                auto fab_arr = mfab1d.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    fab_arr(i,0,0,n) = line[n + ncomp*i];
                });
            }
            break;
        case 1:
            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                auto fab_arr = mfab1d.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    fab_arr(0,j,0,n) = line[n + ncomp*j];
                });
            }
            break;
        case 2:
            for (MFIter mfi(mfab1d, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                auto fab_arr = mfab1d.array(mfi);
                AMREX_FOR_4D(bx, ncomp, i, j, k, n,
                {
                    fab_arr(0,0,k,n) = line[n + ncomp*k];
                });
            }
            break;
    }

    std::string line_plot_file{"line_plot"};

    const std::string& plotfilename = amrex::Concatenate(line_plot_file, m_nstep);

    amrex::Print() << "  Writing line plot file " << plotfilename << " at time " << m_cur_time << std::endl;

    Vector<std::string> pltscaVarsName;

    pltscaVarsName.push_back("<u>");
    pltscaVarsName.push_back("<v>");
    pltscaVarsName.push_back("<w>");
    pltscaVarsName.push_back("<T>");

    pltscaVarsName.push_back("<u'u'>");
    pltscaVarsName.push_back("<u'v'>");
    pltscaVarsName.push_back("<u'w'>");
    pltscaVarsName.push_back("<v'v'>");
    pltscaVarsName.push_back("<v'w'>");
    pltscaVarsName.push_back("<w'w'>");

    pltscaVarsName.push_back("<w'u'u'>");
    pltscaVarsName.push_back("<w'u'v'>");
    pltscaVarsName.push_back("<w'u'w'>");
    pltscaVarsName.push_back("<w'v'v'>");
    pltscaVarsName.push_back("<w'v'w'>");
    pltscaVarsName.push_back("<w'w'w'>");

    pltscaVarsName.push_back("<T'u'>");
    pltscaVarsName.push_back("<T'v'>");
    pltscaVarsName.push_back("<T'w'>");
    pltscaVarsName.push_back("<nu+nu_t>");

    //fixme todo add Rij, qj, nu_SGS https://a2ehfm-ecp.slack.com/archives/C3V26K34G/p1522245519000450

    int level_step = 0;

    // in debug mode amrex output will complain if these are not the same size
    AMREX_ALWAYS_ASSERT(mfab1d.nComp() == pltscaVarsName.size());

    amrex::WriteSingleLevelPlotfile(plotfilename, mfab1d, pltscaVarsName, geom[0], m_cur_time, level_step);

    WriteJobInfo(plotfilename);
    
}
