#ifndef CHANNELFLOWMAP_H
#define CHANNELFLOWMAP_H

#include "amr-wind/core/MeshMap.H"

namespace amr_wind {
namespace channel_map {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class ChannelFlowMap : public MeshMap::Register<ChannelFlowMap>
{
public:
    static const std::string identifier() { return "ChannelFlowMap"; }

    explicit ChannelFlowMap();

    ~ChannelFlowMap() override = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int, const amrex::Geometry&);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int, const amrex::Geometry&);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&);

private:
    //! User input parameters
    amrex::Vector<amrex::Real> m_beta{{0.0, 3.0, 0.0}};

    amrex::Real m_eps{1e-11};
};

} // namespace channel_map
} // namespace amr_wind

#endif /* CHANNELFLOWMAP_H */
