#include <AMReX_Print.H>
#include "amr-wind/utilities/tagging/TerrainElevation.H"
#include "amr-wind/utilities/io_utils.H"

namespace amr_wind {

void TerrainElevation::load_from_file(const std::string& filename)
{
    ioutils::read_flat_grid_file(filename, m_x, m_y, m_z);
    m_loaded = !m_x.empty() && !m_y.empty() && !m_z.empty();
    amrex::Print() << "TerrainElevation loaded from file: " << filename << "\n";
}

std::shared_ptr<TerrainElevation>& TerrainElevation::get_instance()
{
    static std::shared_ptr<TerrainElevation> instance;
    return instance;
}

void TerrainElevation::ensure_loaded(const std::string& filename)
{
    auto& ptr = get_instance();
    if (!ptr) {
        ptr = std::make_shared<TerrainElevation>();
    }
    if (!ptr->is_loaded()) {
        ptr->load_from_file(filename);
    }
}

} // namespace amr_wind
