#ifndef SUBSETSAMPLER_H
#define SUBSETSAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace sampling {

/** Sample a subset of the entire domain
 *  \ingroup sampling
 *
 * Defines a rectangular subset of the entire domain, and outputs the
 * physical coordinates, computational coordinates, and data at each
 * point within the defined subdomain.  The sellected subdomain is
 * defined in the physical coodinate system in the x direction by
 * `xlim`, in y by `ylim`, and z by `zlim`.  The refinement level of
 * the subdomain is given by `level` (optional, default is 0).
 *
 * # Example of how it can be used: 
 * sampling.labels                          = psubset
 * sampling.psubset.type                    = SubsetSampler
 * sampling.psubset.level                   = 0
 * sampling.psubset.xlim                    = 0 0.5
 * sampling.psubset.ylim                    = 0 0.5
 * sampling.psubset.zlim                    = 0 0.5
 *
 */
class SubsetSampler : public SamplerBase::Register<SubsetSampler>
{
public:
    static const std::string identifier() { return "SubsetSampler"; }

    explicit SubsetSampler(const CFDSim&);

    virtual ~SubsetSampler();

    void initialize(const std::string& key) override;

    void sampling_locations(SampleLocType&) const override;

    void define_netcdf_metadata(const ncutils::NCGroup&) const override;

    void populate_netcdf_metadata(const ncutils::NCGroup&) const override;


    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    int id() const override { return m_id; }
    int& id() override { return m_id; }
    int num_points() const override { return m_npts; }

private:
    const amr_wind::CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;

    Field& m_velocity;
    Field& m_mesh_fac_cc;
    Field& m_nu_coord_cc;

    amrex::Vector<vs::Vector> m_pos_nu_vec;
    amrex::Vector<vs::Vector> m_pos_comp_vec;

    int m_level{0};
    amrex::Vector<amrex::Real> m_xlim;
    amrex::Vector<amrex::Real> m_ylim;
    amrex::Vector<amrex::Real> m_zlim;

    amrex::Vector<int> m_proc_offsets;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace sampling
} // namespace amr_wind

#endif /* SUBSETSAMPLER_H */
