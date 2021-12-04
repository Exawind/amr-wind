#ifndef CHANNELFLOWSCALING_H
#define CHANNELFLOWSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Channel flow scaling mesh map
 *  \ingroup mesh_map
 */
class ChannelFlowScaling : public MeshMap::Register<ChannelFlowScaling>
{
public:
    static const std::string identifier() { return "ChannelFlowScaling"; }

    explicit ChannelFlowScaling(const CFDSim& sim);

    virtual ~ChannelFlowScaling() = default;

    //! Construct the mesh scaling field
    void create_map(int, const amrex::Geometry&) override;

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int, const amrex::Geometry&) override{};

private:
    Field& m_mesh_scale_fac_cc;
    Field& m_mesh_scale_fac_nd;

    //! User input parameters
    amrex::Vector<int> n_cell{{0, 0, 0}};
    amrex::Vector<amrex::Real> m_beta{{0.0, 3.0, 0.0}};
    amrex::Vector<int> m_map{{0, 1, 0}};
};

} // namespace amr_wind

#endif /* CHANNELFLOWSCALING_H */
