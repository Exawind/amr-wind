#ifndef SyntheticTurbulence_H
#define SyntheticTurbulence_H

#include <string>
#include <cmath>
#include <memory>

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

namespace synth_turb {
class MeanProfile
{
public:
  MeanProfile(double ref_vel, int shear_dir=2)
    : m_ref_vel(ref_vel),
      m_shear_dir(shear_dir)
  {}

  virtual ~MeanProfile() = default;

  virtual double operator()(double /* height */) const
  { return m_ref_vel; }

  inline double reference_velocity() const { return m_ref_vel; }

  inline int shear_dir() const { return m_shear_dir; }

protected:
  const double m_ref_vel;
  const int m_shear_dir;
};

} // namespace synth_turb

struct SynthTurbTraits
{
  static constexpr int n_dim_max{3};
  using ArrayType = double[n_dim_max];
  using IntArrayType = int[n_dim_max];
  using TransMatType = double[n_dim_max][n_dim_max];
  // using StructField = Kokkos::View<double*, Kokkos::LayoutRight, MemSpace>;
  using StructField = std::vector<double>;
};

struct SynthTurbData
{
  // Dimensions of the box
  vs::VectorT<int> box_dims;

  // Length of the boxes in each direction
  vs::Vector box_len;

  vs::Vector dx;

  // Reference point for the turbulence box. Reference is the mid-point of the
  // turbulence grid at the plane where it is injected into the CFD flow field.
  vs::Vector origin;

  // Transformation matrix to convert from global coordinate system to local
  // coordinate system
  vs::Tensor tr_mat;

  // Perturbation velocities (2, ny, nz)
  std::vector<double> uvel;
  std::vector<double> vvel;
  std::vector<double> wvel;

  // Indices of the two planes stored in the data arrays
  int ileft;
  int iright;
};



/**
 *  \defgroup synth_turb Synthetic Turbulence
 *
 *  SyntheticTurbulence contains functions to inject turbulence into a
 *  CFD simulation. It reads a file to get the turbulent velocity on a
 *  grid and populates a source term for the momentum equation to
 *  achieve the desired turbulence characteristics at a specified
 *  plane in the CFD domain.
 *
 * \ingroup physics
 */

class SyntheticTurbulence: public Physics::Register<SyntheticTurbulence>
{
public:

  static const std::string identifier() { return "SyntheticTurbulence"; }

  struct NCBoxTurb
  {
    std::string filename;

    // NetCDF file ID
    int ncid;

    // Dimensions
    int s_dim, x_dim, y_dim, z_dim;

    // Perturbation velocity field IDs
    int uid, vid, wid;

    // Box length and grid size IDs
    int boxlen_id, dx_id;

    // Scale and divergence correction IDs
    int scale_id, divcorr_id;
  };

  SyntheticTurbulence(const CFDSim& sim);

  SyntheticTurbulence() = delete;
  SyntheticTurbulence(const SyntheticTurbulence&) = delete;
  SyntheticTurbulence& operator=(const SyntheticTurbulence&) = delete;

  virtual ~SyntheticTurbulence() = default;

  void initialize_fields(int level, const amrex::Geometry& geom) override;

  void post_init_actions() override {}

  void post_regrid_actions() override {}

  void pre_advance_work() override;

  void post_advance_work() override {}

private:
  void initialize();

  void update();

  const amr_wind::SimTime& m_time;
  const FieldRepo& m_repo;
  const amrex::AmrCore& m_mesh;
  Field& m_velocity;
  Field& m_density;
  Field& m_turb_force;

  // Turbulence file read information
  NCBoxTurb m_turb_file;

  // Turbulence box data
  SynthTurbData m_turb_grid;

  std::unique_ptr<synth_turb::MeanProfile> m_wind_profile;

  amrex::Real m_grid_spacing;
  amrex::Real m_epsilon;
  amrex::Real m_gauss_scaling;

  amrex::Real m_time_offset{0.0};

  bool m_is_init{true};
};


} // amr_wind

#endif /* SyntheticTurbulence_H */
