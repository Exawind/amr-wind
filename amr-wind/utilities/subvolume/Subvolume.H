#ifndef SUBVOLUME_H
#define SUBVOLUME_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/subvolume/SubvolumeBase.H"
#include <AMReX_PlotFileUtil.H>

/**
 *  \defgroup subvolume Mesh data subvolume utilities
 *  Mesh data subvolume utilities
 *
 *  Subvolume is a data-retrieval tool designed to obtain partitions (chunks)
 * directly from the AMR-Wind computational mesh. It is distinct from the
 * sampling utilities because it does not use particles and is therefore
 * confined to the actual cell locations; however, it is another way to sample
 * data from the simulation. In essence, subvolume outputs are like plotfiles
 * but only contain a portion of the domain.
 *
 *  \ingroup utilities
 */

namespace amr_wind::subvolume {

class Subvolume : public PostProcessBase::Register<Subvolume>
{
public:
    static std::string identifier() { return "Subvolume"; }

    Subvolume(CFDSim& /*sim*/, std::string /*label*/);

    ~Subvolume() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs for chunks, output timing
    void initialize() override;

    //! Do any work that needs to happen after regrid
    void post_regrid_actions() override;

    //! Get chunks and output to disk
    void post_advance_work() override;

    void write_subvolume();

    const amrex::Vector<std::string>& var_names() const { return m_var_names; }

protected:
private:
    CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<SubvolumeBase>> m_subvolumes;

    //! List of variable names for output
    amrex::Vector<std::string> m_var_names;

    //! List of fields to be sampled for this collection of probes
    amrex::Vector<Field*> m_fields;

    //! List of IntFields to be sampled for this collection of probes
    amrex::Vector<IntField*> m_int_fields;

    //! List of derived fields to be sampled for this collection of probes
    std::unique_ptr<DerivedQtyMgr> m_derived_mgr;

    /** Name of this subvolume object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories.
     */
    const std::string m_label;

    // number of field components
    int m_ncomp{0};

    // number of int field components
    int m_nicomp{0};

    // number of derived field components
    int m_ndcomp{0};

    //! Frequency of data output
    int m_out_freq{100};

    //! Delay number of timestep before output
    int m_out_delay{0};
};

} // namespace amr_wind::subvolume

#endif /* SUBVOLUME_H */
