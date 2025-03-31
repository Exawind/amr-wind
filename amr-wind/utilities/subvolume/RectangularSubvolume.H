#ifndef LINESAMPLER_H
#define LINESAMPLER_H

#include "amr-wind/utilities/subvolume/SubvolumeBase.H"

namespace amr_wind::subvolume {

/** Get subvolume defined by 3D rectangle
 *  \ingroup subvolume
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segments as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class RectangularSubvolume
    : public SubvolumeBase::Register<RectangularSubvolume>
{
public:
    static std::string identifier() { return "Rectangular"; }

    explicit RectangularSubvolume(const CFDSim& /*sim*/);

    ~RectangularSubvolume() override;

    /** Read user inputs and initialize the subvolume object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Check if inputs have a valid correspondence in the mesh and save
    //! quantities for output steps
    void evaluate_inputs() override;

    void post_regrid_actions() override { evaluate_inputs(); }

    //! Name of this subvolume object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this subvolume object
    std::string subvolumetype() const override { return identifier(); }

    int lev() const override { return m_lev_for_sub; }

    amrex::BoxArray box_array() const override { return m_ba; }

protected:
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<int> m_npts_vec;
    amrex::Vector<amrex::Real> m_dx_vec;
    amrex::Vector<int> m_chunk_size_vec;

    int m_lev_for_sub{0};

    amrex::BoxArray m_ba;

    std::string m_label;
};

} // namespace amr_wind::subvolume

#endif /* LINESAMPLER_H */
