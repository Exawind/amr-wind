#ifndef CHANNELFLOWMAP_H
#define CHANNELFLOWMAP_H

#include "amr-wind/core/MeshMap.H"

namespace amr_wind::channel_map {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class ChannelFlowMap : public MeshMap::Register<ChannelFlowMap>
{
public:
    static std::string identifier() { return "ChannelFlowMap"; }

    explicit ChannelFlowMap();

    ~ChannelFlowMap() override = default;

    //! Construct the mesh scaling field
    void create_map(int /*lev*/, const amrex::Geometry& /*geom*/) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int /*lev*/, const amrex::Geometry& /*geom*/);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int /*lev*/, const amrex::Geometry& /*geom*/);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int /*lev*/, const amrex::Geometry& /*geom*/);

private:
    //! User input parameters
    amrex::Vector<amrex::Real> m_beta{{0.0, 3.0, 0.0}};

    amrex::Real m_eps{1e-11};
};

} // namespace amr_wind::channel_map

#endif /* CHANNELFLOWMAP_H */
