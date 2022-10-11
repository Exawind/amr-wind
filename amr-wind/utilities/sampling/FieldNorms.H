#ifndef FIELDNORMS_H
#define FIELDNORMS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"

namespace amr_wind::field_norms {

/** field norm object
 *  \ingroup Field Norm
 *
 *  A concrete impelementation of the post-processing interface that deals with
 *  norms of fields. FIXME: need to update this
 *
 */
class FieldNorms : public PostProcessBase::Register<FieldNorms>
{
public:
    static std::string identifier() { return "FieldNorms"; }

    FieldNorms(CFDSim& /*sim*/, std::string /*label*/);

    ~FieldNorms() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    //! Interpolate fields at a given timestep and output to disk
    void post_advance_work() override;

    void post_regrid_actions() override {}

    //! Write sampled data in binary format
    void impl_write_native();

    //! calculate the L2 norm of a given field and component
    static amrex::Real L2_norm(amr_wind::Field& field, const int comp);

    const amrex::Vector<std::string>& var_names() const { return m_var_names; }

private:
    //! prepare ASCII file and directory
    void prepare_ascii_file();

    //! Output sampled data in ASCII format
    void write_ascii();

    //!  process all fields and calculate norms
    void process_field_norms();

    //! Reference to the CFD sim
    CFDSim& m_sim;

    //! List of variable names for output
    amrex::Vector<std::string> m_var_names;

    //! List holding norms for all fields and their components
    amrex::Vector<amrex::Real> m_fnorms;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! filename for ASCII output
    std::string m_out_fname;

    //! Frequency of data sampling and output
    int m_out_freq{10};

    //! width in ASCII output
    int m_width{18};

    //! precision in ASCII output
    int m_precision{12};
};

} // namespace amr_wind::field_norms

#endif /* FIELDNORMS_H */
