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
    int lev, amrex::Real ,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("FieldRepo::make_level_from_coarse");
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> fact(new amrex::FArrayBoxFactory());
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder);

    allocate_field_data(ba, dm, *ldata, *fact);

    // TODO: Fill patch logic

    m_leveldata[lev] = std::move(ldata);
    m_factory[lev] = std::move(fact);

    m_is_initialized = true;
}

void FieldRepo::remake_level(
    int lev, amrex::Real ,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("FieldRepo::remake_level");
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> fact(new amrex::FArrayBoxFactory());
    std::unique_ptr<LevelDataHolder> ldata(new LevelDataHolder);

    allocate_field_data(ba, dm, *ldata, *fact);

    // TODO: Fill patch logic

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
        auto found = m_fields.find(name);
        if (found != m_fields.end()) *m_fields[name];
    }

    // Only allow states to be at times and not half steps
    if (nstates > (static_cast<int>(FieldState::NM1) + 1)) {
        amrex::Abort("Invalid number of states specified for field: " + name);
    }

    // Create the field data structures
    std::shared_ptr<FieldInfo> finfo(
        new FieldInfo(name, ncomp, ngrow, nstates, floc));
    for (int i = 0; i < nstates; ++i) {
        const auto fstate = static_cast<FieldState>(i);
        const std::string fname =
            field_impl::field_name_with_state(name, fstate);
        std::unique_ptr<Field> field(new Field(*this, fname, finfo, fstate));

        finfo->m_states[fstate] = field.get();

        // If declare field is called after mesh has been initialized create field multifabs
        if (m_is_initialized)
            allocate_field_data(*field);

        m_fields[fname] = std::move(field);
    }

    return *m_fields[name];
}

Field& FieldRepo::get_field(
    const std::string& name, const FieldState fstate)
{
    const auto fname = field_impl::field_name_with_state(name, fstate);
    const auto found = m_fields.find(fname);
    if (found == m_fields.end()) {
        amrex::Abort("Cannot find field: " + name);
    }

    return *m_fields[fname];
}

void FieldRepo::allocate_field_data(
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    for (auto& it: m_fields) {
        const auto& fname = it.first;
        const auto& field = *(it.second);

        auto ba1 = amrex::convert(ba, field_impl::index_type(field.field_location()));
        std::unique_ptr<amrex::MultiFab> mfab(new amrex::MultiFab);
        mfab->define(
            ba1, dm, field.num_comp(), field.num_grow(), amrex::MFInfo(), factory);

        auto& fdata = level_data.m_data[fname];
        fdata = std::move(mfab);
    }
}

void FieldRepo::allocate_field_data(
    int lev,
    Field& field,
    LevelDataHolder& level_data,
    const amrex::FabFactory<amrex::FArrayBox>& factory)
{
    std::unique_ptr<amrex::MultiFab> mfab(new amrex::MultiFab);
    auto ba = amrex::convert(
        m_mesh.boxArray(lev), field_impl::index_type(field.field_location()));
    mfab->define(
        ba, m_mesh.DistributionMap(lev), field.num_comp(), field.num_grow(),
        amrex::MFInfo(), factory);

    const std::string& fname = field.name();
    auto& fdata = level_data.m_data[fname];
    fdata = std::move(mfab);
}

void FieldRepo::allocate_field_data(Field& field)
{
    for (int lev=0; lev <= m_mesh.finestLevel(); ++lev) {
        allocate_field_data(lev, field, *m_leveldata[lev], *m_factory[lev]);
    }
}

} // namespace amr_wind
