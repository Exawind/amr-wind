//
//  ThirdMomentAveraging.H
//  amr-wind
//


#ifndef ThirdMomentAveraging_H
#define ThirdMomentAveraging_H


#include <AMReX_AmrCore.H>
#include "DirectionSelector.H"
#include "Field.H"
#include "FieldPlaneAveraging.H"

namespace amr_wind {

class ThirdMomentAveraging
{
public:
    
    ThirdMomentAveraging(FieldPlaneAveraging& pa1,
                         FieldPlaneAveraging& pa2,
                         FieldPlaneAveraging& pa3);

    ~ThirdMomentAveraging(){};

    /** evaluate third moment at specific location for both components */
    amrex::Real line_average_interpolated(amrex::Real x, int comp1, int comp2, int comp3) const;
    /** evaluate third moment at specific location for each component */
    amrex::Real line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate third moment at specific cell for both components */
    amrex::Real line_average_cell(int cell, int comp1, int comp2, int comp3) const;
    /** evaluate third moment at specific cell for each component */
    amrex::Real line_average_cell(int cell, int comp) const;

    /** change precision of text file output */
    void set_precision(int p) {m_precision = p;};

private:

    int m_num_moments; /** outer product of components */
    amrex::Vector<amrex::Real> m_third_moments_line; /** line storage for fluctuations */

    int m_precision = 4; /** precision for line plot text file */

    FieldPlaneAveraging& m_plane_average1;
    FieldPlaneAveraging& m_plane_average2;
    FieldPlaneAveraging& m_plane_average3;

public: // public for GPU

    /** fill line storage with averages */
    template <typename IndexSelector>
    void compute_average(
        const IndexSelector& idxOp,
        const amrex::MultiFab& mfab1,
        const amrex::MultiFab& mfab2,
        const amrex::MultiFab& mfab3);
};

}

#endif /* ThirdMomentAveraging_H */
