#ifndef FieldPlaneAveraging_H
#define FieldPlaneAveraging_H

#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/CFDSim.H"
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
template<typename FType>
class FPlaneAveraging
{
public:

    /**
    *  \param field_in [in] Field to be averaged
    *  \param time [in] Time instance to determine output frequencies
    *  \param axis_in [in] Direction along which planes are computed
    *  \param compute_deriv [in] Should the derivative of the averages be computed
    */
    FPlaneAveraging(
      const FType& field_in, const amr_wind::SimTime& time, int axis_in
      , bool compute_deriv=false);
    
    ~FPlaneAveraging() = default;

    void operator()();

    /** evaluate line average at specific location for any average component */
    amrex::Real line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate line average at specific cell for any average component */
    amrex::Real line_average_cell(int cell, int comp) const;

    /** evaluate line average derivative at specific location for any average
        component */
    amrex::Real line_derivative_interpolated(amrex::Real x, int comp) const;
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
    void line_average(int comp, amrex::Vector<amrex::Real> & l_vec);   
    const amrex::Vector<amrex::Real>& line_centroids()
    {
        return m_line_xcentroid;
    };

    const FType& field() const { return m_field; };

protected:
    int m_ncomp; /** number of average components */

    amrex::Vector<amrex::Real>
        m_line_average; /** line storage for the average velocity and tracer
                           variables */
    //! line storage for the derivative of average field
    amrex::Vector<amrex::Real> m_line_deriv; 
    
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

    const FType& m_field;
    const SimTime& m_time;
    const int m_axis;
    const bool m_comp_deriv;

public: // public for GPU
    /** fill line storage with averages */
    template <typename IndexSelector>
    void
    compute_averages(const IndexSelector& idxOp, const amrex::MultiFab& mfab);
    
    /** fill derivatives of line storage */
    void compute_line_derivatives();
};

using FieldPlaneAveraging = FPlaneAveraging<Field>;

/** Specific application of FieldPlaneAveraging to the velocity field
 *  \ingroup statistics we_abl
 */
class VelPlaneAveraging: public FieldPlaneAveraging
{
public:
    
    VelPlaneAveraging(CFDSim& sim, int axis_in);
    
    ~VelPlaneAveraging() = default;

    void operator()();

private:

    amrex::Vector<amrex::Real>
    m_line_hvelmag_average; /** line storage for the average horizontal velocity magnitude */
    //! line storage for the derivative of average horizontal velocity magnitude
    amrex::Vector<amrex::Real> m_line_hvelmag_deriv;

public: // public for GPU
    /** fill line storage with horizontal velocity magnitude averages */
    template <typename IndexSelector>
    void
    compute_hvelmag_averages(const IndexSelector& idx_op,
                             const int h1_idx,
                             const int h2_idx,                             
                             const amrex::MultiFab& mfab);

    /** return vector containing horizontal velocity magnitude average */
    const amrex::Vector<amrex::Real>& line_hvelmag_average() { return m_line_hvelmag_average; };
    
    /** evaluate line haverage at specific location for horizontal velocity magnitude */
    amrex::Real line_hvelmag_average_interpolated(amrex::Real x) const;
    /** evaluate line average at specific cell for horizontal velocity magnitude */
    amrex::Real line_hvelmag_average_cell(int cell) const;

    /** compute derivatives of horizontal velocity magnitude */
    void compute_line_hvelmag_derivatives();
    /** evaluate line average derivative at specific location for horizontal velocity 
        magnitude */
    amrex::Real line_hvelmag_derivative_interpolated(amrex::Real x) const;
    /** evaluate derivative of a line average at specific cell horizontal velocity 
        magnitude */
    amrex::Real line_hvelmag_derivative_of_average_cell(int ind) const;

    void output_line_average_ascii(
        std::string filename, int step, amrex::Real time);
    void output_line_average_ascii(int step, amrex::Real time);

};
    

} // namespace amr_wind

#endif /* FieldPlaneAveraging_H */
