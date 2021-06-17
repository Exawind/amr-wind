#ifndef CONSTANTSCALING_H
#define CONSTANTSCALING_H

#include "amr-wind/core/MeshMap.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Coonstant scaling mesh map
 *  \ingroup mesh_map
 */
class ConstantScaling : public MeshMap::Register<ConstantScaling>
{
public:
    static const std::string identifier() { return "ConstantScaling"; }

    explicit ConstantScaling(const CFDSim& sim);

    virtual ~ConstantScaling() = default;

    //! Construct the mesh scaling field
    void create_map(int finest_level) override;

private:
    Field& m_mesh_scale_fac;

    //! Factor to scale the mesh by
    amrex::Vector<amrex::Real> m_fac{{0.0, 0.0, 0.0}};
};

} // namespace amr_wind

#endif /* CONSTANTSCALING_H */
