#include "Field.H"
#include "FieldRepo.H"
#include "FieldFillPatchOps.H"

namespace amr_wind {

FieldInfo::FieldInfo(
    const std::string& basename,
    const int ncomp,
    const int ngrow,
    const int nstates,
    const FieldLoc floc)
    : m_basename(basename)
    , m_ncomp(ncomp)
    , m_ngrow(ngrow)
    , m_nstates(nstates)
    , m_floc(floc)
    , m_bc_values(AMREX_SPACEDIM*2, amrex::Vector<amrex::Real>(ncomp, 0.0))
    , m_bc_values_dview(ncomp * AMREX_SPACEDIM * 2)
    , m_bcrec(ncomp)
    , m_bcrec_d(ncomp)
{}

FieldInfo::~FieldInfo() = default;

void FieldInfo::copy_bc_to_device() noexcept
{
    amrex::Vector<amrex::Real> h_data(m_ncomp * AMREX_SPACEDIM * 2);

    // Copy data to a flat array for transfer to device
    {
        amrex::Real* hp = h_data.data();
        for (const auto& v: m_bc_values) {
            for (const auto& x: v)
                *(hp++) = x;
        }
    }

    // Transfer BC values to device
    amrex::Real* ptr = m_bc_values_dview.data();
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(
        ptr, h_data.data(), sizeof(amrex::Real) * h_data.size());
#else
    std::memcpy(ptr, h_data.data(), sizeof(amrex::Real)*h_data.size());
#endif

    for (int i=0; i < AMREX_SPACEDIM*2; ++i) {
        m_bc_values_d[i] = ptr;
        ptr += m_ncomp;
    }

#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(
        m_bcrec_d.data(), m_bcrec.data(),
        sizeof(amrex::BCRec) * AMREX_SPACEDIM);
#else
    std::memcpy(m_bcrec_d.data(), m_bcrec.data(), sizeof(amrex::BCRec) * AMREX_SPACEDIM);
#endif
}

Field::Field(
    FieldRepo& repo,
    const std::string& name,
    const std::shared_ptr<FieldInfo>& info,
    const unsigned fid,
    const FieldState state)
    : m_repo(repo), m_name(name), m_info(info), m_id(fid), m_state(state)
{
}

Field::~Field() = default;

Field& Field::state(const FieldState fstate)
{
    auto& fstates = m_info->m_states;
    const int fint = static_cast<int>(fstate);
    return *fstates[fint];
}

const Field& Field::state(const FieldState fstate) const
{
    auto& fstates = m_info->m_states;
    const int fint = static_cast<int>(fstate);
    return *fstates[fint];
}

amrex::MultiFab& Field::operator()(int lev) noexcept
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    return  m_repo.get_multifab(m_id, lev);
}

const amrex::MultiFab& Field::operator()(int lev) const noexcept
{
    BL_ASSERT(lev < m_repo.num_active_levels());
    return m_repo.get_multifab(m_id, lev);
}

amrex::Vector<amrex::MultiFab*> Field::vec_ptrs() noexcept
{
    const int nlevels = m_repo.num_active_levels();
    amrex::Vector<amrex::MultiFab*> ret;
    ret.reserve(nlevels);
    for (int lev = 0; lev < nlevels; ++lev) {
        ret.push_back(&m_repo.get_multifab(m_id, lev));
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
            &m_repo.get_multifab(m_id, lev)));
    }
    return ret;
}

void Field::fillpatch(int lev, amrex::Real time, amrex::MultiFab& mfab)
{
    BL_ASSERT(m_info->m_fillpatch_op);
    auto& fop = *(m_info->m_fillpatch_op);

    fop.fillpatch(lev, time, mfab);
}

void Field::fillpatch_from_coarse(int lev, amrex::Real time, amrex::MultiFab& mfab)
{
    BL_ASSERT(m_info->m_fillpatch_op);
    auto& fop = *(m_info->m_fillpatch_op);

    fop.fillpatch_coarse(lev, time, mfab);
}

void Field::fillpatch(amrex::Real time)
{
    BL_ASSERT(m_info->m_fillpatch_op);
    auto& fop = *(m_info->m_fillpatch_op);
    const int nlevels = m_repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        fop.fillpatch(lev, time, m_repo.get_multifab(m_id, lev));
    }
}

void Field::fillphysbc(amrex::Real time)
{
    BL_ASSERT(m_info->m_fillpatch_op);
    auto& fop = *(m_info->m_fillpatch_op);
    const int nlevels = m_repo.num_active_levels();
    for (int lev=0; lev < nlevels; ++lev) {
        fop.fillphysbc(lev, time, m_repo.get_multifab(m_id, lev));
    }
}

} // namespace amr_wind
