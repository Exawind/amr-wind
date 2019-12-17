#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#else
#include <AMReX_FArrayBox.H>
#endif

void incflo_predict_vels_on_faces (amrex::Box const& ubx,
                                   amrex::Box const& vbx,
                                   amrex::Box const& wbx,
                                   amrex::Array4<amrex::Real> const& u,
                                   amrex::Array4<amrex::Real> const& v,
                                   amrex::Array4<amrex::Real> const& w,
                                   amrex::Array4<amrex::Real const> const& vcc);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}
