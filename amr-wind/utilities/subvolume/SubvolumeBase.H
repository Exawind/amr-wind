#ifndef SUBVOLUMEBASE_H
#define SUBVOLUMEBASE_H

#include <AMReX_RealVect.H>
#include "amr-wind/core/Factory.H"

namespace amr_wind {

class CFDSim;

namespace subvolume {

/** Abstract representation of subvolume constraints to sample flow data
 *  \ingroup subvolume
 *
 *  This class defines the interface by which the positions of the subvolumes
 *  are defined and access from the Subvolume class. All subvolume definitions
 *  must be sub-classes of this base class.
 */
class SubvolumeBase : public Factory<SubvolumeBase, CFDSim&>
{
public:
    static std::string base_identifier() { return "SubvolumeBase"; }

    ~SubvolumeBase() override = default;

    //! Name used to refer to this sampler (e.g., myline1)
    virtual std::string label() const = 0;
    virtual std::string& label() = 0;

    //! Class name of this sampler (e.g., Rectangular)
    virtual std::string subvolumetype() const = 0;

    //! Read inputs and perform initialization actions
    virtual void initialize(const std::string& key) = 0;

    //! Check inputs and prepare for output
    virtual void evaluate_inputs() = 0;

    //! Run actions after regrid
    virtual void post_regrid_actions() { evaluate_inputs(); }

    //! Level that the subvolume corresponds to
    virtual int lev() const { return 0; }

    //! BoxArray that the subvolume corresponds to
    virtual amrex::BoxArray box_array() const { return amrex::BoxArray(0); }
};

} // namespace subvolume
} // namespace amr_wind

#endif /* SUBVOLUMEBASE_H */
