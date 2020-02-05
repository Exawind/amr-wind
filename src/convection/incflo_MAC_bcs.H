#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <AMReX_FArrayBox.H>

void incflo_set_mac_bcs (amrex::Box const& domain,
                         amrex::Box const& ubx,
                         amrex::Box const& vbx,
                         amrex::Box const& wbx,
                         amrex::Array4<amrex::Real> const& u,
                         amrex::Array4<amrex::Real> const& v,
                         amrex::Array4<amrex::Real> const& w,
                         amrex::Array4<amrex::Real const> const& vel,
                         amrex::Vector<amrex::BCRec> const& h_bcrec,
                         amrex::BCRec const* d_bcrec);

