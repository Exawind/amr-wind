#ifndef ISOPLANESAMPLER_H
#define ISOPLANESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind {
namespace sampling {

/** Sample data on uniform grids
 *  \ingroup isosampling
 *
 *  Defines probe locations on a two-dimensional plane defined by two vectors
 *  `axis1` and `axis2` that originate from the point `origin`. The length of
 *  the vectors determine the extents of the parallelogram on which the probes
 *  are defined. The grid points are uniformly spaced along each axis as
 *  determined by the number of points defined in `num_points` (a list of two
 *  integers, one for each axis).
 *
 *  With no additional inputs, the implementation creates one plane of grid
 *  points such that `origin` lines on the plane defined by the points. The
 *  orientation vector determines the search direction for the isosampling
 *  routine, with each plane sampler having a single orientation vector.
 */
class IsoPlaneSampler : public SamplerBase::Register<IsoPlaneSampler>
{
public:
    static const std::string identifier() { return "IsoPlaneSampler"; }

    explicit IsoPlaneSampler(const CFDSim&);

    virtual ~IsoPlaneSampler();

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType&) const override;
    //! Populate and return a vector of probe orientations to be sampled
    //  (the orientation is constant for a each point in a plane sampler)
    void sampling_orientations(SampleLocType&) const override;

    void define_netcdf_metadata(const ncutils::NCGroup&) const override;
    void populate_netcdf_metadata(const ncutils::NCGroup&) const override;

    //! Name of this sampling object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Unique identifier for this set of probe locations
    int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of probe locations along the line
    int num_points() const override { return m_npts; }

private:
    amrex::Vector<amrex::Real> m_axis1;
    amrex::Vector<amrex::Real> m_axis2;
    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<amrex::Real> m_oris;
    amrex::Vector<int> m_npts_dir;

    std::string m_label;

    int m_id{-1};
    int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* ISOPLANESAMPLER_H */
