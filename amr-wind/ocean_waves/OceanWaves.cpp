#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/OceanWavesModel.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

#include <algorithm>

namespace amr_wind::ocean_waves {

OceanWaves::OceanWaves(CFDSim& sim)
    : m_sim(sim)
    , m_ow_levelset(sim.repo().declare_field("ow_levelset", 1, 3, 1))
    , m_ow_vof(sim.repo().declare_field("ow_vof", 1, 2, 1))
    , m_ow_velocity(
          sim.repo().declare_field("ow_velocity", AMREX_SPACEDIM, 3, 1))
{
    if (!sim.physics_manager().contains("MultiPhase")) {
        amrex::Abort("OceanWaves requires Multiphase physics to be active");
    }
    m_ow_levelset.set_default_fillpatch_bc(sim.time());
    m_ow_vof.set_default_fillpatch_bc(sim.time());
    m_ow_velocity.set_default_fillpatch_bc(sim.time());
}

OceanWaves::~OceanWaves() = default;

void OceanWaves::pre_init_actions()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::pre_init_actions");
    amrex::ParmParse pp(identifier());

    std::string label;
    pp.query("label", label);
    const std::string& tname = label;
    const std::string& prefix = identifier() + "." + tname;
    amrex::ParmParse pp1(prefix);

    std::string type;
    pp.query("type", type);
    pp1.query("type", type);
    AMREX_ALWAYS_ASSERT(!type.empty());

    m_owm = OceanWavesModel::create(type, m_sim, tname, 0);

    const std::string default_prefix = identifier() + "." + type;
    ::amr_wind::utils::MultiParser inp(default_prefix, prefix);

    m_owm->read_inputs(inp);
}

void OceanWaves::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::initialize_fields");
    m_owm->init_waves(level, geom);
}

void OceanWaves::post_init_actions()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::post_init_actions");
    m_owm->update_relax_zones();
}

void OceanWaves::post_regrid_actions() {}

void OceanWaves::pre_advance_work()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::pre_advance_work");
}

void OceanWaves::post_advance_work()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::post_init_actions");
    relaxation_zones();
}

/** Update ocean waves relaxation zones
 *
 */
void OceanWaves::relaxation_zones()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::update_relaxation_zones");
    m_owm->update_relax_zones();
    m_owm->apply_relax_zones();
}

void OceanWaves::prepare_outputs()
{
    const std::string out_dir_prefix = "post_processing/ocean_waves";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }

    m_owm->prepare_outputs(sname);
}

} // namespace amr_wind::ocean_waves