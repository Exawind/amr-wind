#include "amr-wind/utilities/DerivedQtyDefs.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/fvm/fvm.H"
#include "amr-wind/utilities/io_utils.H"

namespace amr_wind {
namespace derived {

VorticityMag::VorticityMag(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.size() == 0u);
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
    AMREX_ALWAYS_ASSERT(args.size() == 0u);
}

void QCriterion::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto q_crit = fld.subview(scomp, 1);
    fvm::q_criterion(q_crit, m_vel);
}

StrainRateMag::StrainRateMag(
    const FieldRepo& repo, const std::vector<std::string>& args)
    : m_vel(repo.get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(args.size() == 0u);
}

void StrainRateMag::operator()(ScratchField& fld, const int scomp)
{
    AMREX_ASSERT(fld.num_comp() > (scomp));
    auto srate = fld.subview(scomp, 1);
    fvm::strainrate(srate, m_vel);
}

Gradient::Gradient(const FieldRepo& repo, const std::vector<std::string>& args)
{
    AMREX_ALWAYS_ASSERT(args.size() == 1u);
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
    AMREX_ALWAYS_ASSERT(args.size() == 1u);
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
    AMREX_ALWAYS_ASSERT(args.size() == 1u);
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
    AMREX_ALWAYS_ASSERT(nargs > 1u);
    m_fld = &repo.get_field(args[0]);
    AMREX_ALWAYS_ASSERT(static_cast<int>(nargs - 1) < m_fld->num_comp());

    m_ncomp = nargs - 1;
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

} // namespace derived
} // namespace amr_wind
