#ifndef PLANESAMPLER_H
#define PLANESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind::sampling {

/** Sample data on uniform grids
 *  \ingroup sampling
 *
 *  Defines probe locations on a two-dimensional plane defined by two vectors
 *  `axis1` and `axis2` that originate from the point `origin`. The length of
 *  the vectors determine the extents of the parallelogram on which the probes
 *  are defined. The grid points are uniformly spaced along each axis as
 *  determined by the number of points defined in `num_points` (a list of two
 *  integers, one for each axis).
 *
 *  With no additional inputs, the implementation creates one plane of grid
 *  points such that `origin` lines on the plane defined by the points. User can
 *  choose to create additional parallel planes by providing an `offset` (a list
 *  of translation distances) from origin along a `normal` vector. Note that if
 *  `offset` is specified, then the implementation will not create a default
 *  plane at `origin`, the user must include a zero translation offset if
 *  sampling on the plane at `origin` is desired.
 */
class PlaneSampler : public SamplerBase::Register<PlaneSampler>
{
public:
    static std::string identifier() { return "PlaneSampler"; }

    explicit PlaneSampler(const CFDSim& /*unused*/);

    ~PlaneSampler() override;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*locs*/) const override;

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;

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
    amrex::Vector<amrex::Real> m_normal{0.0, 0.0, 0.0};
    amrex::Vector<amrex::Real> m_poffsets;
    amrex::Vector<int> m_npts_dir;

    std::string m_label;

    int m_id{-1};
    int m_npts{0};
};

} // namespace amr_wind::sampling

#endif /* PLANESAMPLER_H */
