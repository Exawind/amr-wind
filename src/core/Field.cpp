#include "Field.H"
#include "FieldRepo.H"

namespace amr_wind {

Field::Field(
    FieldRepo& repo,
    const std::string& name,
    const std::shared_ptr<FieldInfo>& info,
    const FieldState state)
    : m_repo(repo), m_name(name), m_info(info), m_state(state)
{}

Field& Field::state(const FieldState fstate)
{
    auto& fstates = m_info->m_states;
    auto found = fstates.find(fstate);
    if (found == fstates.end())
        amrex::Abort("Cannot find requested state for field: " + m_name);
    return *fstates[fstate];
}

const Field& Field::state(const FieldState fstate) const
{
    auto& fstates = m_info->m_states;
    auto found = fstates.find(fstate);
    if (found == fstates.end())
        amrex::Abort("Cannot find requested state for field: " + m_name);
    return *fstates[fstate];
}

amrex::MultiFab& Field::operator()(int lev) noexcept
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    auto* mfab = m_repo.get_multifab(name(), lev);
    BL_ASSERT(mfab != nullptr);
    return *mfab;
}

const amrex::MultiFab& Field::operator()(int lev) const noexcept
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    auto* mfab = m_repo.get_multifab(name(), lev);
    BL_ASSERT(mfab != nullptr);
    return *mfab;
}

amrex::Vector<amrex::MultiFab*> Field::vec_ptrs() noexcept
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::MultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(m_repo.get_multifab(m_name, lev));
    }
    return ret;
}

amrex::Vector<const amrex::MultiFab*> Field::vec_const_ptrs() const noexcept
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<const amrex::MultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(static_cast<const amrex::MultiFab*>(
            m_repo.get_multifab(m_name, lev)));
    }
    return ret;
}

} // namespace amr_wind
