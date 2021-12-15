#ifndef FREESURFACE_H
#define FREESURFACE_H

#include <memory>

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"

/**
 *  \defgroup findinterface Multiphase-sampling utilities
 *  Multiphase-sampling utilities
 *
 *  findinterface uses level_mask to find the location of the free surface
 *  between liquid and gas phases, given a user-defined 2D grid. It supports
 *  output in ascii as well as NetCDF format.
 *
 *  \ingroup utilities
 */

namespace amr_wind {
namespace free_surface {

/** Collection of data sampling objects
 *  \ingroup findinterface
 *
 *  A concrete implementation of the post-processing interface that deals with
 *  data probes. Sampling positions are defined with a plane and a mesh-aligned
 *  orientation.
 */
class FreeSurface : public PostProcessBase::Register<FreeSurface>
{
public:
    static const std::string identifier() { return "FreeSurface"; }

    FreeSurface(CFDSim&, const std::string&);

    ~FreeSurface();

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    //! Interpolate fields at a given timestep and output to disk
    void post_advance_work() override;

    void post_regrid_actions() override {}

    //! Output functions for private variables
    //! Number of points
    int num_gridpoints() const { return m_npts; }
    //! Number of instances
    int num_instances() const { return m_ninst; }
    //! Locations
    amrex::Vector<amrex::Array<amrex::Real, 2>> locations() const
    {
        return m_locs;
    }
    //! Outputs (heights)
    amrex::Vector<amrex::Real> heights() const { return m_out; }

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
    void write_ascii();

private:
    CFDSim& m_sim;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! reference to VOF
    const Field& m_vof;

    //! Format of the data output (ascii, netcdf, etc.)
#ifdef AMR_WIND_USE_NETCDF
    std::string m_out_fmt{"netcdf"};
    std::string m_ncfile_name;
#else
    std::string m_out_fmt{"ascii"};
#endif

    //! Number or points on 2D grid in each direction
    amrex::Vector<int> m_npts_dir;
    int m_npts;
    //! Number of instances (possible sampling points per location, like in the
    //! case of a breaking wave)
    int m_ninst{1};

    //! Parameters to set up plane
    amrex::Vector<amrex::Real> m_start, m_end;
    //! Locations of points in 2D grid
    amrex::Vector<amrex::Array<amrex::Real, 2>> m_locs;
    //! Output coordinate
    amrex::Vector<amrex::Real> m_out;

    //! Orientation of interface search (always in z)
    const int m_orient = 2;
    const amrex::Vector<int> m_griddim = {0, 1};

    //! Frequency of data sampling and output
    int m_out_freq{100};
};

} // namespace free_surface
} // namespace amr_wind

#endif /* FREESURFACE_H */
