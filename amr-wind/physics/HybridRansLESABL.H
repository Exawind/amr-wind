#ifndef HybridRANSLESABL_H
#define HybridRANSLESABL_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {
namespace hybrid_rans_les_abl {

/** Hybrid RANS/LES with ABL
 *  \ingroup physics
 */
class HybridRANSLESABL : public Physics::Register<HybridRANSLESABL>
{
public:
    static const std::string identifier() { return "HybridRANSLESABL"; }

    HybridRANSLESABL(CFDSim& sim);

    virtual ~HybridRANSLESABL() = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const FieldRepo& m_repo;

    Field& m_tke;
    Field& m_sdr;

    //! Turbulence model constant
    amrex::Real m_Ce{0.1};

    //! Compute sdr field from sfs ke field
    void compute_sdr();
};
} // namespace hybrid_rans_les_abl
} // namespace amr_wind

#endif /* HybridRANSLESABL_H */
