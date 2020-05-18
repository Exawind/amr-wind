//
//  FieldPlaneAveraging.H
//  amr-wind
//

#ifndef FieldPlaneAveraging_H
#define FieldPlaneAveraging_H


#include <AMReX_AmrCore.H>
#include "DirectionSelector.H"
#include "Field.H"

class incflo;

class FieldPlaneAveraging
{
public:
    FieldPlaneAveraging(amr_wind::Field &field, int axis_in = 2);

    ~FieldPlaneAveraging(){};

    void operator()();

    amrex::Real eval_line_average(amrex::Real x, int comp) const; /** evaluate line average at specific location for any average component */
    amrex::Real eval_line_average(int cell, int comp) const; /** evaluate line average at specific for any average component */

    /** output a single text file with averages only */
//    void
//    plot_line_average_text(std::string filename, int step, amrex::Real time);
//    /** output a single text file with averages and fluctuations */
//    void plot_line_text(std::string filename, int step, amrex::Real time);
//    /** output a single binary file with averages and fluctuations */
//    void plot_line_binary(std::string filename, int step, amrex::Real time);
//
//    /** output a line plot file based on the plot type */
//    void plot_line(int step, amrex::Real time, int plot_type = 0);

    /** change precision of text file output */
    void set_precision(int p) {precision = p;};

    amrex::Real get_dx(){return dx;};
    amrex::Real get_xlo(){return xlo;};
    int get_axis(){return axis;};
    int get_level(){return level;};
    int get_line_avg_size() {return line_average.size();};
    int get_ncell_plane(){return ncell_plane;};
    int get_ncell_line(){return ncell_line;};
    amrex::Real * get_line_avg_data() {return line_average.data();};
    amr_wind::Field & get_field() {return field;};

private:

    int ncomp; /** number of average components */

    amrex::Vector<amrex::Real> line_average; /** line storage for the average velocity and tracer variables */
    amrex::Vector<amrex::Real> line_xcentroid; /** line storage for centroids of each cell along a line*/

    amrex::Real dx; /** mesh spacing in axis direction*/
    amrex::Real xlo; /** bottom of domain in axis direction */

    int ncell_plane; /** number of cells in plane */
    int ncell_line; /** number of cells along line */

    int precision = 4; /** precision for line plot text file */
    const int level = 0; /** level for plane averaging for now fixed at level=0 */


    amr_wind::Field& field;
    const int axis;

public: // public for GPU

    /** fill line storage with averages */
    template <typename IndexSelector>
    void avg_line(
        const IndexSelector& idxOp,
        const amrex::MultiFab& mfab);
};

#endif /* FieldPlaneAveraging_H */
