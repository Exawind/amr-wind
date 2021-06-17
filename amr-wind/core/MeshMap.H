#ifndef MESHMAP_H
#define MESHMAP_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/CollMgr.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

namespace amr_wind {

class CFDSim;

/**
 *  \defgroup mesh_map Mesh mapping models
 *
 *  AMR-Wind representation of different mesh mapping models
 *
 *  In AMR-Wind, different mesh mappings are implemented using the MeshMap class.
 */

/** Abstract representation of different physics
 *  \ingroup physics
 *
 *  This class defines an abstract API that represents the notion of some
 *  mesh mapping that will be used to scale the mesh. The most common use-case
 *  for this class is to perform RANS simulations.
 */
class MeshMap : public Factory<MeshMap, CFDSim&>
{
public:
    static const std::string base_identifier() { return "MeshMap"; }

    virtual ~MeshMap() = default;

    //! Construct mesh scaling field
    virtual void create_map(int finest_level) = 0;
};

/** A collection of \ref mesh_map instances
 *
 *  \ingroup mesh_map
 */
class MeshMapMgr : public CollMgr<MeshMapMgr, MeshMap>
{
public:
    //! Return a concrete mesh mapping instance
    template <typename T>
    T& get()
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<T&>(operator()(T::identifier()));
    }

    //! Return a concrete mesh mapping instance
    template <typename T>
    const T& get() const
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<const T&>(operator()(T::identifier()));
    }
};

} // namespace amr_wind

#endif /* MESHMAP_H */
