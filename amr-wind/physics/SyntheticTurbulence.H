#ifndef SyntheticTurbulence_H
#define SyntheticTurbulence_H

#include <string>
#include <cmath>
#include <memory>

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

namespace synth_turb {
class MeanProfile
{
public:
  MeanProfile(double refVel, double refHeight)
    : refVel_(refVel),
      refHeight_(refHeight)
  {}

  virtual ~MeanProfile() = default;

  virtual double operator()(double /* height */) const
  { return refVel_; }

  inline double reference_velocity() const { return refVel_; }

  inline double reference_height() const { return refHeight_; }

protected:
  const double refVel_;
  const double refHeight_;
};

} // namespace synth_turb

struct SynthTurbTraits
{
  static constexpr int NDimMax{3};
  using ArrayType = double[NDimMax];
  using IntArrayType = int[NDimMax];
  using TransMatType = double[NDimMax][NDimMax];
  // using StructField = Kokkos::View<double*, Kokkos::LayoutRight, MemSpace>;
  using StructField = std::vector<double>;
};

struct SynthTurbData
{
  // Dimensions of the box
  SynthTurbTraits::IntArrayType boxDims_;

  // Length of the boxes in each direction
  SynthTurbTraits::ArrayType boxLen_;

  SynthTurbTraits::ArrayType dx_;

  // Reference point for the turbulence box. Reference is the mid-point of the
  // turbulence grid at the plane where it is injected into the CFD flow field.
  SynthTurbTraits::ArrayType origin_;

  // Transformation matrix to convert from global coordinate system to local
  // coordinate system
  SynthTurbTraits::TransMatType trMat_;

  // Perturbation velocities (2, ny, nz)
  SynthTurbTraits::StructField uvel_;
  SynthTurbTraits::StructField vvel_;
  SynthTurbTraits::StructField wvel_;
  // SynthTurbTraits::StructField::HostMirror h_uvel_;
  // SynthTurbTraits::StructField::HostMirror h_vvel_;
  // SynthTurbTraits::StructField::HostMirror h_wvel_;

  // Indices of the two planes stored in the data arrays
  int iLeft_;
  int iRight_;
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
    int sDim, xDim, yDim, zDim;

    // Perturbation velocity field IDs
    int uid, vid, wid;

    // Box length and grid size IDs
    int boxLenid, dxid;

    // Scale and divergence correction IDs
    int scaleid, divCorrid;
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
  NCBoxTurb turbFile_;

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
