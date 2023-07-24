#include "CheckpointFileUtil.H"
#include <AMReX_AsyncOut.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace amrex;

namespace {
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace

CheckpointFileDataImpl::CheckpointFileDataImpl(
    std::string const& chkptfile_name, Vector<std::string> addl_field_names)
    : m_chkptfile_name(chkptfile_name)
{
    // Add field names from input
    const int naf = addl_field_names.size();
    for (int n = 0; n < naf; ++n) {
        m_var_names.push_back(addl_field_names[n]);
    }
    // Header
    std::string File(chkptfile_name + "/Header");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::istringstream is(
        std::string(fileCharPtr.dataPtr()), std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file

    // Title line
    std::getline(is, line);

    // Finest level
    is >> m_finest_level;
    GotoNextLine(is);

    // Step count
    is >> m_nstep;
    GotoNextLine(is);

    // Current time
    is >> m_time;
    GotoNextLine(is);

    // m_time.set_restart_time(nstep, cur_time);

    // Time step size
    is >> m_dt_restart;
    GotoNextLine(is);

    is >> m_dt_nm1;
    GotoNextLine(is);

    is >> m_dt_nm2;
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_prob_lo[i++] = std::stod(word);
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_prob_hi[i++] = std::stod(word);
        }
    }

    // Read in box array to get domain box
    m_nlevels = m_finest_level + 1;
    m_ba.resize(m_nlevels);
    m_prob_domain.resize(m_nlevels);
    for (int lev = 0; lev <= m_finest_level; ++lev) {
        // read in level 'lev' BoxArray from Header
        m_ba[lev].readFrom(is);
        GotoNextLine(is);
        // get minimal box from box array for prob domain
        m_prob_domain[lev] = m_ba[lev].minimalBox();
    }

    // m_var_names
    // add to defaults?

    
    AMREX_ASSERT(m_nlevels > 0 && m_nlevels <= 1000);

    m_mf_name.resize(m_nlevels);
    m_vismf.resize(m_nlevels);
    m_ba.resize(m_nlevels);
    m_dmap.resize(m_nlevels);
    m_ngrow.resize(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        int levtmp, ngrids, levsteptmp;
        Real gtime;
        is >> levtmp >> ngrids >> gtime;
        is >> levsteptmp;
        Real glo[3], ghi[3];
        AMREX_ASSERT(ngrids >= 0 && ngrids < std::numeric_limits<int>::max());
        for (int igrid = 0; igrid < ngrids; ++igrid) {
            for (int idim = 0; idim < m_spacedim; ++idim) {
                is >> glo[idim] >> ghi[idim];
            }
        }
        std::string relname;
        is >> relname;
        m_mf_name[ilev] = m_chkptfile_name + "/" + relname;
        if (m_ncomp > 0) {
            m_vismf[ilev] = std::make_unique<VisMF>(m_mf_name[ilev]);
            m_ba[ilev] = m_vismf[ilev]->boxArray();
            m_dmap[ilev].define(m_ba[ilev]);
            m_ngrow[ilev] = m_vismf[ilev]->nGrowVect();
        }
    }
}

void CheckpointFileDataImpl::syncDistributionMap(
    CheckpointFileDataImpl const& src) noexcept
{
    int nlevs_min = std::min(m_nlevels, src.m_nlevels);
    for (int ilev = 0; ilev < nlevs_min; ++ilev) {
        syncDistributionMap(ilev, src);
    }
}

void CheckpointFileDataImpl::syncDistributionMap(
    int level, CheckpointFileDataImpl const& src) noexcept
{
    if (level <= src.finestLevel() &&
        m_dmap[level].size() == src.DistributionMap(level).size()) {
        m_dmap[level] = src.DistributionMap(level);
    }
}

MultiFab CheckpointFileDataImpl::get(int level) noexcept
{
    MultiFab mf(m_ba[level], m_dmap[level], m_ncomp, m_ngrow[level]);
    VisMF::Read(mf, m_mf_name[level]);
    return mf;
}

MultiFab
CheckpointFileDataImpl::get(int level, std::string const& varname) noexcept
{
    MultiFab mf(m_ba[level], m_dmap[level], 1, m_ngrow[level]);
    auto r = std::find(std::begin(m_var_names), std::end(m_var_names), varname);
    if (r == std::end(m_var_names)) {
        amrex::Abort(
            "CheckpointFileDataImpl::get: varname not found " + varname);
    } else {
        int icomp = static_cast<int>(std::distance(std::begin(m_var_names), r));
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            int gid = mfi.index();
            FArrayBox& dstfab = mf[mfi];
            std::unique_ptr<FArrayBox> srcfab(
                m_vismf[level]->readFAB(gid, icomp));
            dstfab.copy<RunOn::Host>(*srcfab);
        }
    }
    return mf;
}