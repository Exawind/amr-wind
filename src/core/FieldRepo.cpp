#include "FieldRepo.H"

namespace amr_wind {

void FieldRepo::make_new_level_from_scratch(
    int lev, amrex::Real /* time */,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("FieldRepo::make_new_level_from_scratch");
    m_factory[lev].reset(new amrex::FArrayBoxFactory());
    m_leveldata[lev].reset(new LevelDataHolder);

    allocate_field_data(ba, dm, *m_leveldata[lev], *m_factory[lev]);

    m_is_initialized = true;
}

void FieldRepo::make_new_level_from_coarse(
    int lev, amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("FieldRepo::make_level_from_coarse");
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> fact(new amrex::FArrayBoxFactory());
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder);

    allocate_field_data(ba, dm, *ldata, *fact);

    for (auto& field: m_field_vec) {
        if (!field->fillpatch_on_regrid()) continue;

        field->fillpatch_from_coarse(lev, time, ldata->m_mfabs[field->id()], 0);
    }

    m_leveldata[lev] = std::move(ldata);
    m_factory[lev] = std::move(fact);

    m_is_initialized = true;
}

void FieldRepo::remake_level(
    int lev, amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("FieldRepo::remake_level");
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> fact(new amrex::FArrayBoxFactory());
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder);

    allocate_field_data(ba, dm, *ldata, *fact);

    for (auto& field: m_field_vec) {
        if (!field->fillpatch_on_regrid()) continue;

        field->fillpatch(lev, time, ldata->m_mfabs[field->id()], 0);
    }

    m_leveldata[lev] = std::move(ldata);
    m_factory[lev] = std::move(fact);

    m_is_initialized = true;
}

void FieldRepo::clear_level(int lev)
{
    BL_PROFILE("FieldRepo::clear_level");
    m_leveldata[lev].reset();
    m_factory[lev].reset();
}

Field& FieldRepo::declare_field(
    const std::string& name,
    const int ncomp,
    const int ngrow,
    const int nstates,
    const FieldLoc floc)
{
    // If the field is already registered check and return the fields
    {
        auto found = m_fid_map.find(name);
        if (found != m_fid_map.end()) {
            auto& field = *m_field_vec[found->second];

            if ((ncomp != field.num_comp()) ||
                (field.num_time_states() != nstates) ||
                (floc != field.field_location())) {
                amrex::Abort("Attempt to reregister field with inconsistent parameters: "
                             + name);
            }
            return field;
        }
    }

    // Only allow states to be at times and not half steps
    if (nstates > (static_cast<int>(FieldState::NM1) + 1)) {
        amrex::Abort("Invalid number of states specified for field: " + name);
    }

    if (!field_impl::is_valid_field_name(name)) {
        amrex::Abort("Attempt to use reserved field name: " + name);
    }

    // Create the field data structures
    std::shared_ptr<FieldInfo> finfo(
        new FieldInfo(name, ncomp, ngrow, nstates, floc));
    for (int i = 0; i < nstates; ++i) {
        const auto fstate = static_cast<FieldState>(i);
        const std::string fname =
            field_impl::field_name_with_state(name, fstate);
        const unsigned fid = m_field_vec.size();

        // Create new field instance
        std::unique_ptr<Field> field(new Field(*this, fname, finfo, fid, fstate));
        // If declare field is called after mesh has been initialized create field multifabs
        if (m_is_initialized)
            allocate_field_data(*field);

        // Add reference to states lookup
        finfo->m_states[i] = field.get();
        // Store the field instance
        m_field_vec.emplace_back(std::move(field));
        // Name to ID lookup map
        m_fid_map[fname] = fid;
    }

    // We want the New state returned when we have multiple states involved
    return *m_field_vec[m_fid_map[name]];
}

Field& FieldRepo::get_field(
    const std::string& name, const FieldState fstate) const
{
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_fid_map.find(fname);
    if (found == m_fid_map.end()) {
        amrex::Abort("Cannot find field: " + name);
    }

    AMREX_ASSERT(found->second < static_cast<unsigned>(m_field_vec.size()));
    return *m_field_vec[found->second];
}

bool FieldRepo::field_exists(
    const std::string& name, const FieldState fstate) const
{
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_fid_map.find(fname);
    return (found != m_fid_map.end());
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field(
    const std::string& name, const int ncomp, const int nghost, const FieldLoc floc) const
{
    if (!m_is_initialized) {
        amrex::Abort("Scratch field creation is not permitted before mesh is initialized");
    }

    std::unique_ptr<ScratchField> field(new ScratchField(*this, name, ncomp, nghost, floc));

    for (int lev=0; lev <= m_mesh.finestLevel(); ++lev) {
        const auto ba = amrex::convert(
            m_mesh.boxArray(lev), field_impl::index_type(floc));

        field->m_data.emplace_back(
            ba, m_mesh.DistributionMap(lev), ncomp, nghost, amrex::MFInfo(), *m_factory[lev]);
    }
    return field;
}

std::unique_ptr<ScratchField> FieldRepo::create_scratch_field(
    const int ncomp, const int nghost, const FieldLoc floc) const
{
    return create_scratch_field("scratch_field", ncomp, nghost, floc);
}

void FieldRepo::advance_states() noexcept
{
    for (auto& it: m_field_vec) {
        if (it->field_state() != FieldState::New) continue;
        it->advance_states();
    }
}

void FieldRepo::allocate_field_data(
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    auto& mfab_vec = level_data.m_mfabs;

    for (auto& field : m_field_vec) {
        auto ba1 =
            amrex::convert(ba, field_impl::index_type(field->field_location()));

        mfab_vec.emplace_back(
            ba1, dm, field->num_comp(), field->num_grow(), amrex::MFInfo(),
            factory);

        mfab_vec.back().setVal(0.0);
    }
}

void FieldRepo::allocate_field_data(
    int lev,
    Field& field,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    auto& mfab_vec = level_data.m_mfabs;
    AMREX_ASSERT(mfab_vec.size() == field.id());
    const auto ba = amrex::convert(
        m_mesh.boxArray(lev), field_impl::index_type(field.field_location()));

    mfab_vec.emplace_back(
        ba, m_mesh.DistributionMap(lev), field.num_comp(), field.num_grow(),
        amrex::MFInfo(), factory);

    mfab_vec.back().setVal(0.0);
}

void FieldRepo::allocate_field_data(Field& field)
{
    for (int lev=0; lev <= m_mesh.finestLevel(); ++lev) {
        allocate_field_data(lev, field, *m_leveldata[lev], *m_factory[lev]);
    }
}

Field& FieldRepo::create_state(
    Field& infield, const FieldState fstate) noexcept
{
    AMREX_ASSERT((fstate == FieldState::NPH));
    AMREX_ASSERT(!field_exists(infield.base_name(), fstate));

    auto& finfo = infield.m_info;
    const int i = static_cast<int>(fstate);
    const std::string fname = field_impl::field_name_with_state(
        infield.base_name(), fstate);
    const unsigned fid = m_field_vec.size();

    // Create the half state
    std::unique_ptr<Field> field(new Field(*this, fname, finfo, fid, fstate));
    if (m_is_initialized)
        allocate_field_data(*field);

    // Add reference to states lookup
    finfo->m_states[i] = field.get();
    m_field_vec.emplace_back(std::move(field));
    m_fid_map[fname] = fid;

    return *m_field_vec.back();
}

} // namespace amr_wind
