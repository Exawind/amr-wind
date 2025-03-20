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
    Field m_bodyz0;
    Field m_body_height;

    std::string m_body_file;
    std::string m_roughness_file;

    bool m_body_is_waves{false};
    Field const* m_wave_volume_fraction{nullptr};
    Field const* m_wave_negative_elevation{nullptr};
    std::string m_wave_volume_fraction_name{"vof"};
    std::string m_wave_negative_elevation_name{"wave_neg_elev"};
};

} // namespace amr_wind::superellipsebody

#endif
