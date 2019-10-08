#ifndef _INCFLO_DIVOP_CONV_HPP_
#define _INCFLO_DIVOP_CONV_HPP_

#include <incflo.H>


void 
incflo_apply_eb_redistribution ( Box& bx,
                                 MultiFab& conv,
                                 MultiFab& divc,
                                 MFIter* mfi,
                                 const int conv_comp,
                                 const int ncomp,
                                 const EBCellFlagFab& flags_fab,
                                 const MultiFab* volfrac,
                                 Box& domain,
                                 const int cyclic_x,
                                 const int cyclic_y,
                                 const int cyclic_z,
                                 const Real* dx);
#endif
