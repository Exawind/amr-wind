#ifndef AMREX_SUPER_ELLIPSE_BODY_H
#define AMREX_SUPER_ELLIPSE_BODY_H

#include "amr-wind/physics/PhysicsBase.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::superellipsebody {

class SuperEllipseBody : public PhysicsBase
{
  public:
    SuperEllipseBody(CFDSim& sim);
    ~SuperEllipseBody() override = default;

    std::string identifier() const override { return "SuperEllipseBody"; }

    void initialize_fields(int level, const amrex::Geometry& geom) override;
    void post_init_actions() override;
    void pre_advance_work() override;
    void post_regrid_actions() override;

  private:
    void convert_waves_to_body_fields();

    CFDSim& m_sim;
    FieldRepo& m_repo;
    Mesh& m_mesh;

    amrex::IntMultiFab m_body_blank;
    amrex::IntMultiFab m_body_drag;

    std::string m_body_file;

};

} // namespace amr_wind::superellipsebody

#endif
