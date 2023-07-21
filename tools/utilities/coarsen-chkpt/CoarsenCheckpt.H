#ifndef COARSENCHECKPT_H
#define COARSENCHECKPT_H

#include "amr-wind/incflo.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind {
namespace tools {

class IOManager_Mod : public IOManager
{
public:
    explicit IOManager_Mod(CFDSim& /*sim*/);

    void read_checkpoint_fields_offset(
        const std::string& restart_file,
        const amrex::Vector<amrex::BoxArray>& ba_chk,
        const amrex::Vector<amrex::DistributionMapping>& dm_chk,
        const amrex::IntVect& rep,
        const int off);

private:
    CFDSim& m_sim;

    //! Final list of fields for restart
    amrex::Vector<Field*> m_chk_fields;

    //! Flag indicating whether we should allow missing restart fields
    bool m_allow_missing_restart_fields{false};
};

class CoarsenCheckpt : public incflo
{
public:
    CoarsenCheckpt();

    virtual ~CoarsenCheckpt() = default;

    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

    void read_chkpt_add_baselevel();

    void coarsen_chkpt_file();

    void run_utility();

    IOManager_Mod& io_manager() { return *m_io_mgr; }
    const IOManager_Mod& io_manager() const { return *m_io_mgr; }

private:
    amrex::Vector<amrex::BoxArray> m_orig_ba;

    std::unique_ptr<IOManager_Mod> m_io_mgr;
};

} // namespace tools
} // namespace amr_wind

#endif /* COARSENCHECKPT_H */
