#ifndef COARSENCHECKPT_H
#define COARSENCHECKPT_H

#include "amr-wind/incflo.H"

namespace amr_wind {
namespace tools {

class CoarsenCheckpt : public incflo
{
public:
    CoarsenCheckpt();

    virtual ~CoarsenCheckpt() = default;

    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

    void read_chkpt_file();

    void coarsen_chkpt_file();

    void run_utility();

private:
    amrex::Vector<amrex::BoxArray> m_orig_ba;
};

} // namespace tools
} // namespace amr_wind

#endif /* COARSENCHECKPT_H */
