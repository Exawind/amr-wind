#ifndef FieldPlaneAveraging_H
#define FieldPlaneAveraging_H

#include <AMReX_AmrCore.H>
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/SimTime.H"

/**
 *  \defgroup statistics Field statistics
 *  Field statistics
 *
 *  This group contains utilities for performing turbulence averaging and
 *  outputting statistics during wind simulations.
 *
 *  \ingroup utilities
 */

namespace amr_wind {

/** Output average of a field on planes normal to a given direction
 *  \ingroup statistics we_abl
 *
 *  The user can choose a direction (x, y, or z), the default value is the
 *  z-direction. The field is then averaged on planes at the cell-centers at
 *  level 0 (coarsest level) along the specified direction.
 */
class FieldPlaneAveraging
{
public:
    FieldPlaneAveraging(
        amr_wind::Field& field, const amr_wind::SimTime& time, int axis_in = 2);

    ~FieldPlaneAveraging(){};

    void operator()();

    /** evaluate line average at specific location for any average component */
    amrex::Real line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate line average at specific cell for any average component */
    amrex::Real line_average_cell(int cell, int comp) const;

    /** evaluate derivative of a line average at specific cell for any component */
    amrex::Real line_derivative_of_average_cell(int ind, int comp) const;

    void
    output_line_average_ascii(std::string filename, int step, amrex::Real time);
    void output_line_average_ascii(int step, amrex::Real time);

    /** change precision of text file output */
    void set_precision(int p) { m_precision = p; };

    amrex::Real dx() const { return m_dx; };
    amrex::Real xlo() const { return m_xlo; };

    int axis() const { return m_axis; };
    int level() const { return m_level; };
    int ncomp() const { return m_ncomp; };
    int ncell_plane() const { return m_ncell_plane; };
    int ncell_line() const { return m_ncell_line; };
    int last_updated_index() const { return m_last_updated_index; };

    const amrex::Vector<amrex::Real>& line_average() { return m_line_average; };
    const amrex::Vector<amrex::Real>& line_centroids()
    {
        return m_line_xcentroid;
    };

    amr_wind::Field& field() const { return m_field; };

private:
    int m_ncomp; /** number of average components */

    amrex::Vector<amrex::Real>
        m_line_average; /** line storage for the average velocity and tracer
                           variables */
    amrex::Vector<amrex::Real> m_line_xcentroid; /** line storage for centroids
                                                    of each cell along a line*/

    amrex::Real m_dx;  /** mesh spacing in axis direction*/
    amrex::Real m_xlo; /** bottom of domain in axis direction */

    int m_ncell_plane; /** number of cells in plane */
    int m_ncell_line;  /** number of cells along line */

    int m_precision = 4; /** precision for line plot text file */
    const int m_level =
        0; /** level for plane averaging for now fixed at level=0 */
    int m_last_updated_index = -1; /** keep track of the last time index that
                                      the operator was called */

    amr_wind::Field& m_field;
    const SimTime& m_time;
    const int m_axis;

public: // public for GPU
    /** fill line storage with averages */
    template <typename IndexSelector>
    void
    compute_averages(const IndexSelector& idxOp, const amrex::MultiFab& mfab);
};

} // namespace amr_wind

#endif /* FieldPlaneAveraging_H */
