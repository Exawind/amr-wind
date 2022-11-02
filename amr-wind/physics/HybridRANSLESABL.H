#ifndef HybridRANSLESABL_H
#define HybridRANSLESABL_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::hybrid_rans_les_abl {

/** Hybrid RANS/LES with ABL
 *  \ingroup physics
 */
class HybridRANSLESABL : public Physics::Register<HybridRANSLESABL>
{
public:
    static std::string identifier() { return "HybridRANSLESABL"; }

    explicit HybridRANSLESABL(const CFDSim& sim);

    ~HybridRANSLESABL() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    //! Compute sdr field from sfs ke field
    void compute_sdr_impl();

private:
    const CFDSim& m_sim;

    Field* m_tke{nullptr};
    Field* m_sdr{nullptr};

    //! Turbulence model constant
    amrex::Real m_Ce{0.1};
};
} // namespace amr_wind::hybrid_rans_les_abl

#endif /* HybridRANSLESABL_H */
