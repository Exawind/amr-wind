//
//  Godunov.H
//  amr-wind
//

//

#ifndef Godunov_H
#define Godunov_H

#include "FieldRepo.H"


namespace godunov {

    void compute_convective_term (amr_wind::FieldRepo& repo,
                                  amrex::Real dt,
                                  const amr_wind::FieldState fstate,
                                  bool constant_density,
                                  bool advect_tracer,
                                  bool godunov_ppm);

    void compute_advection(int lev, amrex::Box const& bx, int ncomp,
                           amrex::Array4<amrex::Real> const& dqdt,
                           amrex::Array4<amrex::Real const> const& q,
                           amrex::Array4<amrex::Real const> const& umac,
                           amrex::Array4<amrex::Real const> const& vmac,
                           amrex::Array4<amrex::Real const> const& wmac,
                           amrex::Array4<amrex::Real const> const& fq,
                           amrex::BCRec const* pbc, int const* iconserv,
                           amrex::Real* p,
                           amrex::Vector<amrex::Geometry> geom,
                           amrex::Real dt,
                           bool godunov_ppm);



    
}


#endif /* Godunov_H */
