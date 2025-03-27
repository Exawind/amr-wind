#ifndef SuperEllipseBody_H
#define SuperEllipseBody_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::superellipsebody {

class SuperEllipseBody : public Physics::Register<SuperEllipseBody>
{
  public:
    static std::string identifier() { return "SuperEllipseBody"; }

    explicit SuperEllipseBody(CFDSim& sim);

    ~SuperEllipseBody() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}
    
    void post_init_actions() override {};

    void post_regrid_actions() override;

    void pre_advance_work() override {};

    void post_advance_work() override {}

  private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Blanking field for SuperEllipseBody
    IntField&  m_body_blank;
    IntField&  m_terrain_drag;

    vs::Vector m_loc{0.0, 0.0, 0.0};

    vs::Vector m_dim{1.0, 1.0, 1.0};
    
    vs::Tensor m_orient = vs::Tensor::identity();

    std::string m_body_file;

};

} // namespace amr_wind::superellipsebody

#endif
