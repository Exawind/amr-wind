#ifndef WAVEENERGY_H
#define WAVEENERGY_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"

namespace amr_wind {
namespace wave_energy {

/** wave energy object
 *  \ingroup Wave Energy
 *
 *  A concrete implementation of the post-processing interface that deals with
 *  integrating total wave energy. This routine calculates the kinetic energy
 *  (KE) and potential energy (PE) according to the following definitions. It
 *  allows for a user-defined potential energy offset, as well (PE_0), intended
 *  to give the result PE = 0 for an undisturbed interface.
 *
 *  KE = 0.5 * (Integral over liquid phase) \rho * u^2 dV
 *  PE = (Integral over liquid phase) \rho * g * z_liq dV + PE_0
 *
 *  Note: the output quantities are not normalized.
 */
class WaveEnergy : public PostProcessBase::Register<WaveEnergy>
{
public:
    static std::string identifier() { return "WaveEnergy"; }

    WaveEnergy(CFDSim& /*sim*/, std::string /*label*/);

    ~WaveEnergy() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and get information needed for calculations
    void initialize() override;

    //! Integrate energy components and output to file
    void post_advance_work() override;

    void post_regrid_actions() override {}

    //! Calculate the sum of stated energy in liquid phase
    amrex::Real calculate_kinetic_energy();
    amrex::Real calculate_potential_energy();

    //! Output private variables that store energy measurements
    void wave_energy(amrex::Real& ke, amrex::Real& pe) const
    {
        ke = m_wave_kinetic_energy;
        pe = m_wave_potential_energy;
    };

private:
    //! prepare ASCII file and directory
    void prepare_ascii_file();

    //! Output sampled data in ASCII format
    virtual void write_ascii();

    //! store the total wave energy
    amrex::Real m_wave_kinetic_energy{0.0};
    amrex::Real m_wave_potential_energy{0.0};

    //! Reference to the CFD sim
    CFDSim& m_sim;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! reference to velocity
    const Field& m_velocity;
    //! reference to vof
    const Field& m_vof;

    //! gravity vector
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    //! offset for potential energy calculation
    amrex::Real m_pe_off = 0.0;

    //! density of liquid phase
    amrex::Real m_rho1 = 10.0;

    //! filename for ASCII output
    std::string m_out_fname;

    //! Frequency of data sampling and output
    int m_out_freq{10};

    //! width in ASCII output
    int m_width{22};

    //! precision in ASCII output
    int m_precision{12};
};

} // namespace wave_energy
} // namespace amr_wind

#endif /* WAVEENERGY_H */
