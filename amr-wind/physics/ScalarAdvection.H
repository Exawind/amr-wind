#ifndef ScalarAdvection_H
#define ScalarAdvection_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

/** Boussinesq Bubble physics
 *  \ingroup physics
 *
 *  \sa ScalarAdvectionFieldInit
 */
class ScalarAdvection : public Physics::Register<ScalarAdvection>
{
public:
    static std::string identifier() { return "ScalarAdvection"; }

    explicit ScalarAdvection(CFDSim& sim);

    ~ScalarAdvection() override = default;

    //! Initialize the temperature and velocity fields for ScalarAdvection
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

    template <typename Shape>
    void initialize_scalar(const Shape& scalar_function);

private:
    const CFDSim& m_sim;
    const FieldRepo& m_repo;
    Field& m_velocity;
    Field& m_density;

    //! Scalar field
    Field* m_scalar{nullptr};

    //! Initial wave packet location
    amrex::Real m_center{0.1};

    //! Initial wave amplitude
    amrex::Real m_amplitude{1.0};

    //! Initial wave width (standard deviation for Gaussian and Yalla2021)
    amrex::Real m_width{0.01};

    //! Eta parameter for Yalla2021
    amrex::Real m_eta{16 * utils::pi()};
    
    //! Initial wave packet shape
    std::string m_shape{"gaussian"};

    //! Density
    amrex::Real m_rho{1.0};
};

struct Gaussian
{
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
      const amrex::Real /*x*/,
      const amrex::Real /*center*/,
      const amrex::Real /*amplitude*/,
      const amrex::Real /*width*/,
      const amrex::Real /*unused*/) const;
};

struct TopHat
{
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
      const amrex::Real /*x*/,
      const amrex::Real /*center*/,
      const amrex::Real /*amplitude*/,
      const amrex::Real /*width*/,
      const amrex::Real /*unused*/) const;
};

struct Yalla2021
{
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
      const amrex::Real /*x*/,
      const amrex::Real /*center*/,
      const amrex::Real /*amplitude*/,
      const amrex::Real /*width*/,
      const amrex::Real /*eta*/) const;
};

} // namespace amr_wind

#endif /* ScalarAdvection_H */
