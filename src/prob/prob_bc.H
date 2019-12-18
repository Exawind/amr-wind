#ifndef PROB_BC_H_
#define PROB_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct IncfloVelFill
{
    int probtype;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloVelFill (int a_probtype,
                             amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,
                                             AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();
        const BCRec& bc = bcr[bcomp];

        if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::x,Orientation::low)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::x,Orientation::low)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::x,Orientation::low)][2];
        }
        else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::x,Orientation::high)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::x,Orientation::high)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::x,Orientation::high)][2];
        }

        if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::y,Orientation::low)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::y,Orientation::low)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::y,Orientation::low)][2];
            if (32 == probtype)
            {
                Real z = (k+0.5)*(1.0/domain_box.length(2));
                vel(i,j,k,1) *= 6. * z * (1.-z);
            }
        }
        else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::y,Orientation::high)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::y,Orientation::high)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::y,Orientation::high)][2];
        }

        if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::z,Orientation::low)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::z,Orientation::low)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::z,Orientation::low)][2];
        }
        else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
        {
            vel(i,j,k,0) = bcv[Orientation(Direction::z,Orientation::high)][0];
            vel(i,j,k,1) = bcv[Orientation(Direction::z,Orientation::high)][1];
            vel(i,j,k,2) = bcv[Orientation(Direction::z,Orientation::high)][2];
        }
    }
};



#endif
