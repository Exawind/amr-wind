#ifndef PDEBASE_H
#define PDEBASE_H

#include <string>

#include "Factory.H"
#include "incflo_enums.H"
#include "FieldDescTypes.H"
#include "FieldRepo.H"
#include "PDEHelpers.H"
#include "SimTime.H"
#include "CollMgr.H"

namespace amr_wind {

class CFDSim;

namespace pde {

class PDEBase : public Factory<PDEBase, CFDSim&, const int>
{
public:
    virtual ~PDEBase() = default;

    virtual PDEFields& fields() = 0;

    virtual void initialize() = 0;

    virtual void post_regrid_actions() = 0;

    virtual void compute_source_term(const FieldState fstate) = 0;

    virtual void compute_nueff(const FieldState fstate) = 0;

    virtual void compute_diffusion_term(const FieldState fstate) = 0;

    virtual void compute_advection_term(const FieldState fstate) = 0;

    virtual void compute_predictor_rhs(const DiffusionType difftype) = 0;

    virtual void compute_corrector_rhs(const DiffusionType difftype) = 0;

    virtual void solve(const amrex::Real dt) = 0;

    static std::string base_identifier() { return "PDESystem"; }
};

class PDEMgr : public CollMgr<PDEMgr, PDEBase>
{
public:
    PDEMgr(CFDSim& sim);

    ~PDEMgr() = default;

    PDEBase& icns() { return *m_icns; }

    PDEBase& register_icns();

private:
    CFDSim& m_sim;

    std::unique_ptr<PDEBase> m_icns;

    std::string m_scheme;

    bool m_use_godunov{true};
};

} // namespace pde
} // namespace amr_wind

#endif /* PDEBASE_H */
