#include "IntField.H"
#include "FieldRepo.H"

namespace amr_wind {

IntField::IntField(
    FieldRepo& repo,
    const std::string& name,
    const unsigned fid,
    const int ncomp,
    const int ngrow,
    const FieldLoc floc)
    : m_repo(repo)
    , m_name(name)
    , m_id(fid)
    , m_ncomp(ncomp)
    , m_ngrow(ngrow)
    , m_floc(floc)
{}

amrex::iMultiFab& IntField::operator()(int lev) noexcept
{
    AMREX_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_int_fab(m_id, lev);
}

const amrex::iMultiFab& IntField::operator()(int lev) const noexcept
{
    AMREX_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_int_fab(m_id, lev);
}

amrex::Vector<amrex::iMultiFab*> IntField::vec_ptrs() noexcept
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::iMultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(&m_repo.get_int_fab(m_id, lev));
    }
    return ret;
}

amrex::Vector<const amrex::iMultiFab*> IntField::vec_const_ptrs() const noexcept
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<const amrex::iMultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(static_cast<const amrex::iMultiFab*>(
                          &m_repo.get_int_fab(m_id, lev)));
    }
    return ret;
}

} // namespace amr_wind
