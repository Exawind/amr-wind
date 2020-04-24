#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include <incflo.H>
#include <Physics.H>
#include "console_io.H"

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
        HeaderFile << m_time.time_index() << "\n";
        HeaderFile << m_time.current_time() << "\n";
        HeaderFile << m_time.deltaT() << "\n";
        HeaderFile << m_time.deltaTNm1() << "\n";
        HeaderFile << m_time.deltaTNm2() << "\n";

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
    BL_PROFILE("amr-wind::incflo::WriteCheckPointFile()")

    const std::string& checkpointname = amrex::Concatenate(m_check_file, m_time.time_index());

    amrex::Print() << "Writing checkpoint " << checkpointname << std::endl;

    amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

    bool is_checkpoint = true;
    WriteHeader(checkpointname, is_checkpoint);
    WriteJobInfo(checkpointname);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Write(velocity()(lev),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "velocity"));

        VisMF::Write(density()(lev),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "density"));

        VisMF::Write(temperature()(lev),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "tracer"));

        VisMF::Write(grad_p()(lev),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gradp"));

        VisMF::Write(pressure()(lev),
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p"));
    }
}

void incflo::ReadCheckpointFile()
{
    BL_PROFILE("amr-wind::incflo::ReadCheckpointFile()")

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

    int nstep;
    amrex::Real cur_time;
    // Step count
    is >> nstep;
    GotoNextLine(is);

    // Current time
    is >> cur_time;
    GotoNextLine(is);

    m_time.set_restart_time(nstep, cur_time);

    // Time step size
    is >> m_time.deltaT();
    GotoNextLine(is);

    is >> m_time.deltaTNm1();
    GotoNextLine(is);

    is >> m_time.deltaTNm2();
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

        MakeNewLevelFromScratch(lev, m_time.current_time(), ba, dm);
    }

    /***************************************************************************
     * Load fluid data                                                         *
     ***************************************************************************/

    // Load the field data
    for(int lev = 0; lev <= finest_level; ++lev)
    {
        VisMF::Read(velocity()(lev),
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "velocity"));

        VisMF::Read(density()(lev),
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "density"));

        VisMF::Read(temperature()(lev),
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "tracer"));


        VisMF::Read(grad_p()(lev),
                    amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "gradp"));

        VisMF::Read(pressure()(lev),
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

        FullPathJobInfoFile += "/amr_wind_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        amr_wind::io::print_banner(jobInfoFile);

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
    BL_PROFILE("amr-wind::incflo::WritePlotFile()")

    if (m_plt_vort or m_plt_divu or m_plt_forcing) {
        IntVect ng(1);
        velocity().fillpatch(m_time.new_time(), ng);
        density().fillpatch(m_time.new_time(), ng);
        temperature().fillpatch(m_time.new_time(), ng);
    }

    const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_time.time_index());

    amrex::Print() << "Writing plotfile " << plotfilename << " at time " << m_time.current_time() << std::endl;

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
    if (m_plt_tracer) ncomp += 1;

    // Pressure
    if(m_plt_p) ++ncomp;

    // Apparent viscosity
    if(m_plt_eta) ++ncomp;

    // Vorticity
    if(m_plt_vort) ++ncomp;

    // Forcing terms in velocity update
    if(m_plt_forcing) ncomp += 3;

    // Magnitude of the rate-of-strain tensor 
    if(m_plt_strainrate) ++ncomp;

    // Magnitude of the stress tensor 
    if(m_plt_stress) ++ncomp;

    // Divergence of velocity field
    if(m_plt_divu) ++ncomp;

    Vector<MultiFab> mf(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
    }

    Vector<std::string> pltscaVarsName;
    int icomp = 0;
    if (m_plt_velx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], velocity()(lev), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velx");
        ++icomp;
    }
    if (m_plt_vely) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], velocity()(lev), 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("vely");
        ++icomp;
    }
    if (m_plt_velz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], velocity()(lev), 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("velz");
        ++icomp;
    }
    if (m_plt_gpx) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], grad_p()(lev), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpx");
        ++icomp;
    }
    if (m_plt_gpy) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], grad_p()(lev), 1, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpy");
        ++icomp;
    }
    if (m_plt_gpz) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], grad_p()(lev), 2, icomp, 1, 0);
        }
        pltscaVarsName.push_back("gpz");
        ++icomp;
    }
    if (m_plt_rho) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], density()(lev), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("density");
        ++icomp;
    }
    if (m_plt_tracer) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], temperature()(lev), 0, icomp, 1, 0);
        }
        pltscaVarsName.push_back("tracer0");
        ++icomp;
    }
    if (m_plt_p) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::average_node_to_cellcenter(mf[lev], icomp, pressure()(lev), 0, 1);
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
            ComputeVorticity(lev, m_time.current_time(), vort, velocity()(lev));
        }
        pltscaVarsName.push_back("vort");
        ++icomp;
    }
    if (m_plt_forcing) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(mf[lev], icns().fields().src_term(lev), 0, icomp, AMREX_SPACEDIM, 0);
        }
        pltscaVarsName.push_back("forcing_x");
        pltscaVarsName.push_back("forcing_y");
        pltscaVarsName.push_back("forcing_z");
        icomp += 3;
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

    AMREX_ALWAYS_ASSERT(ncomp == static_cast<int>(pltscaVarsName.size()));

    // This needs to be defined in order to use amrex::WriteMultiLevelPlotfile, 
    // but will never change unless we use subcycling. 
    // If we do use subcycling, this should be a incflo class member. 
    Vector<int> istep(finest_level + 1, m_time.time_index());

    // Write the plotfile
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf), 
                                   pltscaVarsName, Geom(), m_time.current_time(), istep, refRatio());
    WriteJobInfo(plotfilename);
}
