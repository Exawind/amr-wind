#ifndef ABLISTATS_H
#define ABLISTATS_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"

namespace amr_wind {

/**
 *  \defgroup abl_istats ABL Integrated Statistics
 *  Data-sampling utilities
 *
 *  ABLiStats contains functions to compute and integrated statistics for ABL
 *  simulations. It supports output in ascii format as well as NetCDF format.
 * 
 * \ingroup we_abl
 */
class ABLiStats
{
public:
    static const std::string identifier() { return "ABLiStats"; }

    ABLiStats(CFDSim&, const ABLWallFunction&);

    virtual ~ABLiStats();

    //! Read user inputs and create the necessary files
    void initialize();

    //! Process fields given timestep and output to disk
    void post_advance_work();

protected:
    //! Output data based on user-defined format
    virtual void process_output();

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

    //! Compute height of capping inversion
    void compute_zi();

private:
    CFDSim& m_sim;
    const ABLWallFunction& m_abl_wall_func;

    //! Format of the data output (ascii, netcdf, etc.)
#ifdef AMR_WIND_USE_NETCDF
    std::string m_out_fmt{"netcdf"};
    std::string m_ncfile_name;
#else
    std::string m_out_fmt{"native"};
#endif

    //! Frequency of data sampling and output
    int m_out_freq{100};

    //! Variable to store capping inversion height
    double m_zi{0.0};
};


}

#endif /* ABLISTATS_H */
