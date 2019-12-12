#include <incflo.H>

using namespace amrex;

void incflo::set_background_pressure (int lev)
{
    auto& ld = *m_leveldata[lev];

    if (probtype == 11) {
        gp0[0] = gp0[1] = gp0[2] = 0.0;
        ld.p0.setVal(0.0);
        use_boussinesq = true;
    } else {
        
        for (MFIter mfi(ld.p0); mfi.isValid(); ++mfi)
        {
            Box const& gbx = mfi.fabbox();
        }
    }
}
