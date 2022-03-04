#include "amr-wind/utilities/DerivedQtyDefs.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/fvm/fvm.H"
#include "amr-wind/utilities/io_utils.H"
#include "AMReX_MFIter.H"

namespace amr_wind {
namespace derived {

VorticityMag::VorticityMag(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.empty());
}

void VorticityMag::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto vort_mag = fld.subview(scomp, 1);
    fvm::vorticity_mag(vort_mag, m_vel);
}

QCriterion::QCriterion(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.empty());
}

void QCriterion::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto q_crit = fld.subview(scomp, 1);
    fvm::q_criterion(q_crit, m_vel);
}

QCriterionNondim::QCriterionNondim(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.empty());
}

void QCriterionNondim::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto q_crit_nd = fld.subview(scomp, 1);
    fvm::q_criterion(q_crit_nd, m_vel, true);
}

StrainRateMag::StrainRateMag(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.empty());
}

void StrainRateMag::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto srate = fld.subview(scomp, 1);
    fvm::strainrate(srate, m_vel);
}

Gradient::Gradient(const FieldRepo& repo, const std::vector<std::string>& args)
{
    AMREX_ALWAYS_ASSERT(args.size() == 1U);
    m_phi = &repo.get_field(args[0]);
}

void Gradient::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
    auto gradphi = fld.subview(scomp, num_comp());
    fvm::gradient(gradphi, *m_phi);
}

Divergence::Divergence(
    const FieldRepo& repo, const std::vector<std::string>& args)
{
    AMREX_ALWAYS_ASSERT(args.size() == 1U);
    m_phi = &repo.get_field(args[0]);
    AMREX_ALWAYS_ASSERT(m_phi->num_comp() == AMREX_SPACEDIM);
}

void Divergence::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
    auto divphi = fld.subview(scomp, num_comp());
    fvm::divergence(divphi, *m_phi);
}

Laplacian::Laplacian(
    const FieldRepo& repo, const std::vector<std::string>& args)
{
    AMREX_ALWAYS_ASSERT(args.size() == 1U);
    m_phi = &repo.get_field(args[0]);
    AMREX_ALWAYS_ASSERT(m_phi->num_comp() == AMREX_SPACEDIM);
}

void Laplacian::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
    auto lapphi = fld.subview(scomp, num_comp());
    fvm::laplacian(lapphi, *m_phi);
}

FieldComponents::FieldComponents(
    const FieldRepo& repo, const std::vector<std::string>& args)
{
    const size_t nargs = args.size();
    AMREX_ALWAYS_ASSERT(nargs > 1U);
    m_fld = &repo.get_field(args[0]);
    AMREX_ALWAYS_ASSERT(static_cast<int>(nargs - 1) < m_fld->num_comp());

    m_ncomp = static_cast<int>(nargs) - 1;
    m_comp.resize(nargs - 1);
    for (size_t i = 1; i < nargs; ++i) {
        m_comp[i - 1] = std::stoi(args[i]);
        AMREX_ALWAYS_ASSERT(
            (m_comp[i - 1] >= 0) && (m_comp[i - 1] < m_fld->num_comp()));
    }
}

void FieldComponents::var_names(amrex::Vector<std::string>& plt_var_names)
{
    amrex::Vector<std::string> names;
    ioutils::add_var_names(names, m_fld->name(), m_fld->num_comp());
    for (auto ic : m_comp) {
        plt_var_names.push_back(names[ic]);
    }
}

void FieldComponents::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
    int dst_comp = scomp;
    for (auto icomp : m_comp) {
        field_ops::copy(fld, *m_fld, icomp, dst_comp, 1, 0);
        ++dst_comp;
    }
}

Multiply::Multiply(
    const FieldRepo& repo, const std::vector<std::string>& args)
{
    const size_t nargs = args.size();
    AMREX_ALWAYS_ASSERT(nargs > 1U);
    AMREX_ALWAYS_ASSERT(nargs == 2);
    m_fld1 = &repo.get_field(args[0]);
    m_fld2 = &repo.get_field(args[1]);

    m_ncomp1 = m_fld1->num_comp();
    m_ncomp2 = m_fld2->num_comp();
    int ncomp = m_ncomp1 * m_ncomp2;
    m_comp.resize(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        m_comp[i] = i;
    }
}

void Multiply::var_names(amrex::Vector<std::string>& plt_var_names)
{
    amrex::Vector<std::string> names1;
    amrex::Vector<std::string> names2;
    ioutils::add_var_names(names1, m_fld1->name(), m_fld1->num_comp());
    ioutils::add_var_names(names2, m_fld2->name(), m_fld2->num_comp());
    for (int i = 0; i < m_ncomp1; ++i) {
        for (int j = 0; j < m_ncomp2; ++j) {
            plt_var_names.push_back(names1[i] + "T" + names2[j]);
        }
    }
}

void Multiply::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() >= (scomp + num_comp()));
    // Initialize destination component index
    int dst_comp = scomp;
    // Step through components of first field
    for (int nn = 0; nn < m_ncomp1; ++nn) {

              const int nlevels = fld.repo().num_active_levels();
              for (int lev = 0; lev < nlevels; ++lev) {
      #ifdef _OPENMP
      #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
      #endif
                  for (amrex::MFIter mfi(fld(lev), amrex::TilingIfNotGPU());
                       mfi.isValid(); ++mfi) {
                      const auto& bx = mfi.tilebox();
                      const auto& f1 = (*m_fld1)(lev).array(mfi);
                      const auto& f2 = (*m_fld2)(lev).array(mfi);
                      const auto& dst_f = fld(lev).array(mfi);

                      // Step through components of second field
                      amrex::ParallelFor(
                          bx, m_ncomp2,
                          [=] AMREX_GPU_DEVICE(
                              int i, int j, int k, int n) noexcept {
                              dst_f(i, j, k, dst_comp + n) =
                                  f1(i, j, k, nn) * f2(i, j, k, n);
                          });
                  }
              }
              // Increment by number of second field components
              dst_comp += m_ncomp2;
        }

}

} // namespace derived
} // namespace amr_wind
