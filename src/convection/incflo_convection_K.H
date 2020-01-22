#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

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
amrex::Real incflo_xslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    if (edlo and i == domlo) {
        dc = (vcc(i+1,j,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i-1,j,k,n))/3.0;
    } else if (edhi and i == domhi) {
        dc = (4.0*vcc(i+1,j,k,n)-3.0*vcc(i,j,k,n)-vcc(i-1,j,k,n))/3.0;
    }
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
amrex::Real incflo_yslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    if (edlo and j == domlo) {
        dc = (vcc(i,j+1,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j-1,k,n))/3.0;
    } else if (edhi and j == domhi) {
        dc = (4.0*vcc(i,j+1,k,n)-3.0*vcc(i,j,k,n)-vcc(i,j-1,k,n))/3.0;
    }
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    if (edlo and k == domlo) {
        dc = (vcc(i,j,k+1,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j,k-1,n))/3.0;
    } else if (edhi and k == domhi) {
        dc = (4.0*vcc(i,j,k+1,n)-3.0*vcc(i,j,k,n)-vcc(i,j,k-1,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}


#ifdef AMREX_USE_EB

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope_eb (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& vcc,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,k).isConnected(-1,0,0) and flag(i,j,k).isConnected(1,0,0)) {
        return incflo_xslope(i,j,k,n,vcc);
    } else {
        return 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope_extdir_eb (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& vcc,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (flag(i,j,k).isConnected(-1,0,0) and flag(i,j,k).isConnected(1,0,0)) {
        return incflo_xslope_extdir(i,j,k,n,vcc,edlo,edhi,domlo,domhi);
    } else {
        return 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope_eb (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& vcc,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,k).isConnected(0,-1,0) and flag(i,j,k).isConnected(0,1,0)) {
        return incflo_yslope(i,j,k,n,vcc);
    } else {
        return 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope_extdir_eb (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& vcc,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (flag(i,j,k).isConnected(0,-1,0) and flag(i,j,k).isConnected(0,1,0)) {
        return incflo_yslope_extdir(i,j,k,n,vcc,edlo,edhi,domlo,domhi);
    } else {
        return 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope_eb (int i, int j, int k, int n,
                              amrex::Array4<amrex::Real const> const& vcc,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,k).isConnected(0,0,-1) and flag(i,j,k).isConnected(0,0,1)) {
        return incflo_zslope(i,j,k,n,vcc);
    } else {
        return 0.0;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope_extdir_eb (int i, int j, int k, int n,
                                     amrex::Array4<amrex::Real const> const& vcc,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    if (flag(i,j,k).isConnected(0,0,-1) and flag(i,j,k).isConnected(0,0,1)) {
        return incflo_zslope_extdir(i,j,k,n,vcc,edlo,edhi,domlo,domhi);
    } else {
        return 0.0;
    }
}

#endif

}
