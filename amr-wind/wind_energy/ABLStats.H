#ifndef ABLISTATS_H
#define ABLISTATS_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/SecondMomentAveraging.H"
#include "amr-wind/utilities/ThirdMomentAveraging.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"

namespace amr_wind {

/**
 *  \defgroup abl_istats ABL Statistics
 *
 *  ABLStats contains functions to compute and statistics for ABL
 *  simulations. It supports output in ascii format as well as NetCDF format.
 * 
 * \ingroup we_abl
 */
class ABLStats
{
public:
    static const std::string identifier() { return "ABLStats"; }

    ABLStats(CFDSim&, const ABLWallFunction&);

    virtual ~ABLStats();

    //! Read user inputs and create the necessary files
    void initialize();

    //! Calculate plane average profiles
    void calc_averages();
    
    //! Process fields given timestep and output to disk
    void post_advance_work();

    //! Compute height of capping inversion
    template <typename h1_dir, typename h2_dir>
    void compute_zi(
        const h1_dir& h1Sel,
        const h2_dir& h2Sel);

    //! Return vel plane averaging instance
    const VelPlaneAveraging& vel_plane_averaging() { return m_pa_vel;} ;

    //! Return instance that handles temperature statistics
    const FieldPlaneAveraging& temperature_plane_stats() const
    { return m_pa_temp; }

protected:
    //! Output data based on user-defined format
    virtual void process_output();

    //! Prepare ASCII file
    virtual void prepare_ascii_file();
    //! Prepare NetCDF metadata
    virtual void prepare_netcdf_file();

    //! Write sampled data into a NetCDF file
    void write_netcdf();

    /** Output sampled data in ASCII format
     *
     *  Note that this should be used for debugging only and not in production
     *  runs as it can have significant impacts on code performance.
     */
    virtual void write_ascii();


private:
    CFDSim& m_sim;
    const ABLWallFunction& m_abl_wall_func;
    Field& m_temperature;
    Field& m_mueff;

    VelPlaneAveraging m_pa_vel;
    FieldPlaneAveraging m_pa_temp;
    FieldPlaneAveraging m_pa_mueff;    
    SecondMomentAveraging m_pa_tu;
    SecondMomentAveraging m_pa_uu;
    ThirdMomentAveraging m_pa_uuu;
    
    //! Format of the data output (ascii, netcdf, etc.)
#ifdef AMR_WIND_USE_NETCDF
    std::string m_out_fmt{"netcdf"};
    std::string m_ncfile_name;
#else
    std::string m_out_fmt{"ascii"};
#endif
    std::string m_ascii_file_name;

    //! Frequency of data sampling and output
    int m_out_freq{100};

    //! Acceleration due to gravity magnitude
    amrex::Real m_gravity{9.81};

    //! Von-Karman constant
    amrex::Real m_kappa{0.41};

    //! Reference surface temperature
    amrex::Real m_ref_theta{300.0};    

    //! Variable to store capping inversion height
    double m_zi{0.0};

    //! Wall-normal direction axis
    int m_normal_dir{2};
    
    //! Cell spacing at the coarsest level
    double m_dn{0.0};
    
    //! Number of cells in the horizontal direction
    size_t m_ncells_h1{0};
    size_t m_ncells_h2{0};
};


}

#endif /* ABLISTATS_H */
