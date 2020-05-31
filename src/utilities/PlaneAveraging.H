//
//  PlaneAveraging.H
//  amr-wind
//


#ifndef PlaneAveraging_H
#define PlaneAveraging_H

#include <AMReX_AmrCore.H>
#include "amr-wind/utilities/DirectionSelector.H"

class incflo;

class PlaneAveraging
{
public:
    PlaneAveraging(int axis_in = 2);

    ~PlaneAveraging(){};

    void operator()(const amrex::Vector<amrex::Geometry>& Geom,
               amrex::Vector<amrex::MultiFab*> const& density,
               amrex::Vector<amrex::MultiFab*> const& velocity,
               amrex::Vector<amrex::MultiFab*> const& mueff,
               amrex::Vector<amrex::MultiFab*> const& tracer);

    /** evaluate line average at specific location for velocity xdir*/
    amrex::Real line_velocity_xdir(amrex::Real x) const
    {
        return eval_line_average(x, u_avg);
    };

    /** evaluate line average at specific location for velocity ydir*/
    amrex::Real line_velocity_ydir(amrex::Real x) const
    {
        return eval_line_average(x, v_avg);
    };

    /** evaluate line average at specific location for velocity zdir*/
    amrex::Real line_velocity_zdir(amrex::Real x) const
    {
        return eval_line_average(x, w_avg);
    };

    /** output a single text file with averages only */
    void
    plot_line_average_text(std::string filename, int step, amrex::Real time);
    /** output a single text file with averages and fluctuations */
    void plot_line_text(std::string filename, int step, amrex::Real time);
    /** output a single binary file with averages and fluctuations */
    void plot_line_binary(std::string filename, int step, amrex::Real time);

    /** output a line plot file based on the plot type */
    void plot_line(int step, amrex::Real time, int plot_type = 0);

    /** change precision of text file output */
    void set_precision(int p)
    {
        precision = p;
    };

private:
    enum average_var {
        u_avg,
        v_avg,
        w_avg,
        T_avg,
        nu_avg,
        last_avg_var = nu_avg
    };
    enum fluctuation_var {
        uu,
        uv,
        uw,
        vv,
        vw,
        ww,
        wuu,
        wuv,
        wuw,
        wvv,
        wvw,
        www,
        Tu,
        Tv,
        Tw,
        last_fluc_var = Tw
    };
    const static int navg =
        last_avg_var + 1; /** number of average components */
    const static int nfluc = last_fluc_var+1; /** number of fluctuation components */

    amrex::Vector<amrex::Real> line_average; /** line storage for the average velocity and tracer variables */
    amrex::Vector<amrex::Real> line_fluctuation; /** line storage for fluctuations */
    amrex::Vector<amrex::Real> line_xcentroid; /** line storage for centroids of each cell along a line*/

    amrex::Real dx; /** mesh spacing in axis direction*/
    amrex::Real xlo; /** bottom of domain in axis direction */
    
    int ncell_plane; /** number of cells in plane */
    int ncell_line; /** number of cells along line */
    
    int precision = 4; /** precision for line plot text file */
    int level = 0; /** level for plane averaging for now fixed at level=0 */


    const int axis;

    amrex::Real eval_line_average(amrex::Real x, int comp) const; /** evaluate line average at specific location for any average component */


public: // public for GPU 


    /** fill line storage with averages and fluctuations */
    template <typename IndexSelector>
    void fill_line(
        const IndexSelector& idxOp,
        const amrex::MultiFab& density,
        const amrex::MultiFab& velocity,
        const amrex::MultiFab& mueff,
        const amrex::MultiFab& tracer);
};

#endif /* PlaneAveraging_H */


