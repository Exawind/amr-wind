//
//  SecondMoment.H
//  amr-wind
//


#ifndef SecondMoment_H
#define SecondMoment_H


#include <AMReX_AmrCore.H>
#include "DirectionSelector.H"
#include "Field.H"
#include "FieldPlaneAveraging.H"

class incflo;

class SecondMoment
{
public:
    
    SecondMoment(FieldPlaneAveraging& pa1, FieldPlaneAveraging& pa2);

    ~SecondMoment(){};

    amrex::Real eval_second_moment(amrex::Real x, int comp1, int comp2) const; /** evaluate second moment at specific location for both components */
    amrex::Real eval_second_moment(amrex::Real x, int comp) const; /** evaluate second moment at specific location for each component */

    /** change precision of text file output */
    void set_precision(int p) {precision = p;};

private:

    int nfluc; /** outer product of components */
    amrex::Vector<amrex::Real> second_moments_line; /** line storage for fluctuations */

    int ncell_plane; /** number of cells in plane */
    int ncell_line;  /** number of cells along line */

    int precision = 4; /** precision for line plot text file */

    int ncomp1;
    int ncomp2;

    FieldPlaneAveraging& pa1;
    FieldPlaneAveraging& pa2;

    const int axis;
    const int level;


public: // public for GPU

    /** fill line storage with averages */
    template <typename IndexSelector>
    void fill_line(
        const IndexSelector& idxOp,
        const amrex::MultiFab& mfab1,
        const amrex::MultiFab& mfab2);
};


#endif /* SecondMoment_H */
