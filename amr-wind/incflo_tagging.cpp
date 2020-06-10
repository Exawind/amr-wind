#include "amr-wind/incflo.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

using namespace amrex;

// tag cells for refinement
// overrides the pure virtual function in AmrCore
void incflo::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    BL_PROFILE("amr-wind::incflo::ErrorEst()");

    static bool first = true;
    static Vector<Real> rhoerr_v, gradrhoerr_v;

    if (first) {
        first = false;
        ParmParse pp("incflo");

        pp.queryarr("rhoerr", rhoerr_v);
        if (rhoerr_v.size() > 0) {
            Real last = rhoerr_v.back();
            rhoerr_v.resize(max_level+1, last);
        }

        pp.queryarr("gradrhoerr", gradrhoerr_v);
        if (gradrhoerr_v.size() > 0) {
            Real last = gradrhoerr_v.back();
            gradrhoerr_v.resize(max_level+1, last);
        }
    }

    const auto   tagval = TagBox::SET;
//    const auto clearval = TagBox::CLEAR;


//    const auto prob_lo = geom[lev].ProbLoArray();

    bool tag_rho = lev < rhoerr_v.size();
    bool tag_gradrho = lev < gradrhoerr_v.size();

    if (tag_gradrho) {
        density().fillpatch(lev,time,density()(lev),1);
    }

    const auto& den = density()(lev);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(den,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto const& tag = tags.array(mfi);

        if (tag_rho or tag_gradrho) 
        {
            Array4<Real const> const& rho = den.const_array(mfi);
            Real rhoerr = tag_rho ? rhoerr_v[lev]: std::numeric_limits<Real>::max();
            Real gradrhoerr = tag_gradrho ? gradrhoerr_v[lev] : std::numeric_limits<Real>::max();
            amrex::ParallelFor(bx,
            [tag_rho,tag_gradrho,rhoerr,gradrhoerr,rho,tag]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (tag_rho and rho(i,j,k) > rhoerr) {
                    tag(i,j,k) = tagval;
                }
                if (tag_gradrho) {
                    Real ax = amrex::Math::abs(rho(i+1,j,k) - rho(i,j,k));
                    Real ay = amrex::Math::abs(rho(i,j+1,k) - rho(i,j,k));
                    Real az = amrex::Math::abs(rho(i,j,k+1) - rho(i,j,k));
                    ax = amrex::max(ax,amrex::Math::abs(rho(i,j,k) - rho(i-1,j,k)));
                    ay = amrex::max(ay,amrex::Math::abs(rho(i,j,k) - rho(i,j-1,k)));
                    az = amrex::max(az,amrex::Math::abs(rho(i,j,k) - rho(i,j,k-1)));
                    if (amrex::max(ax,ay,az) >= gradrhoerr) {
                        tag(i,j,k) = tagval;
                    }
                }
            });
        } 
    }

    for (auto& rc: m_refine_criteria) {
        (*rc)(lev, tags, time, ngrow);
    }

}
