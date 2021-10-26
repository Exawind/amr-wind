#ifndef ISOLINESAMPLER_H
#define ISOLINESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind {
namespace sampling {

/** Sample data along a line
 *  \ingroup isosampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segements as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f] . The
 *  line itself is just the initial location for the isosampling routine, which
 *  searches along a single orientation vector for the points belonging to the
 *  line sampler.
 */
class IsoLineSampler : public SamplerBase::Register<IsoLineSampler>
{
public:
    static const std::string identifier() { return "IsoLineSampler"; }

    explicit IsoLineSampler(const CFDSim&);

    virtual ~IsoLineSampler();

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType&) const override;
    //! Populate and return a vector of probe orientations to be sampled
    //  (the orientation is constant for a each point in a line sampler)
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
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;
    amrex::Vector<amrex::Real> m_oris;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* ISOLINESAMPLER_H */
