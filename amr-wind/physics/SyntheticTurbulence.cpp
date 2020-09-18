#include "amr-wind/physics/SyntheticTurbulence.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

#include "netcdf.h"

namespace amr_wind {
namespace synth_turb {

class LinearShearProfile : public MeanProfile
{
public:
  LinearShearProfile(amrex::Real refVel, amrex::Real refHeight, amrex::Real slope, amrex::Real height)
    : MeanProfile(refVel, refHeight),
      slope_(slope),
      halfHeight_(0.5 * height)
  {}

  virtual ~LinearShearProfile() = default;

  virtual amrex::Real operator()(amrex::Real ht) const
  {
    const amrex::Real relHt = ht - refHeight_;
    if (relHt < -halfHeight_)
      return refVel_ * (1.0 - slope_ * halfHeight_);
    else if (relHt > halfHeight_)
      return refVel_ * (1.0 + slope_ * halfHeight_);
    else
      return refVel_ * (1.0 + slope_ * relHt);
  }

private:
  amrex::Real slope_;
  amrex::Real halfHeight_;
};

class PowerLawProfile : public MeanProfile
{
public:
  PowerLawProfile(amrex::Real refVel, amrex::Real refHeight, amrex::Real alpha)
    : MeanProfile(refVel, refHeight),
      alpha_(alpha)
  {}

  virtual ~PowerLawProfile() = default;

  virtual amrex::Real operator()(amrex::Real height) const override
  {
    return refVel_ * std::pow((height / refHeight_), alpha_);
  }

private:
  const amrex::Real alpha_;
};

} // namespace synth_turb

namespace {
/** Check NetCDF errors and throw runtime errors with a message
 */
inline void check_nc_error(int ierr)
{
  if (ierr != NC_NOERR)
    throw std::runtime_error(
      "SyntheticTurbulence NetCDF Error: " + std::string(nc_strerror(ierr)));
}

/** Parse the NetCDF turbulence database and determine details of the turbulence box.
 *
 *. Initializes the dimensions and grid length, sizes in SynthTurbData. Also
 *  allocates the necessary memory for the perturbation velocities.
 *
 *. @param turbFile Information regarding NetCDF data identifiers
 *. @param turbGrid Turbulence data
 */
void process_nc_file(
  SyntheticTurbulence::NCBoxTurb& turbFile,
  SynthTurbData& turbGrid)
{
  check_nc_error(nc_open(turbFile.filename.c_str(), NC_NOWRITE, &turbFile.ncid));

  size_t ndim, nx, ny, nz;
  check_nc_error(nc_inq_dimid(turbFile.ncid, "ndim", &turbFile.sDim));
  check_nc_error(nc_inq_dimlen(turbFile.ncid, turbFile.sDim, &ndim));
  AMREX_ASSERT(ndim == SynthTurbTraits::NDimMax);

  // Grid dimensions
  check_nc_error(nc_inq_dimid(turbFile.ncid, "nx", &turbFile.xDim));
  check_nc_error(nc_inq_dimlen(turbFile.ncid, turbFile.xDim, &nx));
  check_nc_error(nc_inq_dimid(turbFile.ncid, "ny", &turbFile.yDim));
  check_nc_error(nc_inq_dimlen(turbFile.ncid, turbFile.yDim, &ny));
  check_nc_error(nc_inq_dimid(turbFile.ncid, "nz", &turbFile.zDim));
  check_nc_error(nc_inq_dimlen(turbFile.ncid, turbFile.zDim, &nz));

  turbGrid.boxDims_[0] = nx;
  turbGrid.boxDims_[1] = ny;
  turbGrid.boxDims_[2] = nz;

  // Box lengths and resolution
  check_nc_error(nc_inq_varid(turbFile.ncid, "box_lengths", &turbFile.boxLenid));
  check_nc_error(nc_get_var_double(turbFile.ncid, turbFile.boxLenid, turbGrid.boxLen_));
  check_nc_error(nc_inq_varid(turbFile.ncid, "dx", &turbFile.dxid));
  check_nc_error(nc_get_var_double(turbFile.ncid, turbFile.dxid, turbGrid.dx_));

  // Perturbation velocity info
  check_nc_error(nc_inq_varid(turbFile.ncid, "uvel", &turbFile.uid));
  check_nc_error(nc_inq_varid(turbFile.ncid, "vvel", &turbFile.vid));
  check_nc_error(nc_inq_varid(turbFile.ncid, "wvel", &turbFile.wid));
  nc_close(turbFile.ncid);

  // Create data structures to store the perturbation velocities for two planes
  // [t, t+dt] such that the time of interest is within this interval.
  // turbGrid.uvel_ = SynthTurbTraits::StructField("SynthTurbData::uvel", 2*ny*nz);
  // turbGrid.vvel_ = SynthTurbTraits::StructField("SynthTurbData::vvel", 2*ny*nz);
  // turbGrid.wvel_ = SynthTurbTraits::StructField("SynthTurbData::wvel", 2*ny*nz);
  // turbGrid.h_uvel_ = Kokkos::create_mirror_view(turbGrid.uvel_);
  // turbGrid.h_vvel_ = Kokkos::create_mirror_view(turbGrid.vvel_);
  // turbGrid.h_wvel_ = Kokkos::create_mirror_view(turbGrid.wvel_);
  const size_t gridSize = 2 * ny * nz;
  turbGrid.uvel_.resize(gridSize);
  turbGrid.vvel_.resize(gridSize);
  turbGrid.wvel_.resize(gridSize);
}

/** Load two planes of data that bound the current timestep
 *
 *  The data for the y and z directions are loaded for the entire grid at the two planes
 */
void load_turb_plane_data(
  SyntheticTurbulence::NCBoxTurb& turbFile,
  SynthTurbData& turbGrid,
  const int il, const int ir)
{
  check_nc_error(nc_open(turbFile.filename.c_str(), NC_NOWRITE, &turbFile.ncid));

  size_t start[SynthTurbTraits::NDimMax]{static_cast<size_t>(il), 0, 0};
  size_t count[SynthTurbTraits::NDimMax]{
    2, static_cast<size_t>(turbGrid.boxDims_[1]),
    static_cast<size_t>(turbGrid.boxDims_[2])};
  if ((ir - il) == 1) {
    // two consequtive planes load them in one shot
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.uid, start, count, &turbGrid.uvel_[0]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.vid, start, count, &turbGrid.vvel_[0]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.wid, start, count, &turbGrid.wvel_[0]));
  } else {
    // Load the planes separately
    count[0] = 1;
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.uid, start, count, &turbGrid.uvel_[0]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.vid, start, count, &turbGrid.vvel_[0]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.wid, start, count, &turbGrid.wvel_[0]));

    start[0] = static_cast<size_t>(ir);
    const size_t offset = turbGrid.boxDims_[1] * turbGrid.boxDims_[2];
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.uid, start, count, &turbGrid.uvel_[offset]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.vid, start, count, &turbGrid.vvel_[offset]));
    check_nc_error(nc_get_vara_double(
      turbFile.ncid, turbFile.wid, start, count, &turbGrid.wvel_[offset]));
  }

  // Update left and right indices for future checks
  turbGrid.iLeft_ = il;
  turbGrid.iRight_ = ir;

  nc_close(turbFile.ncid);
}

/** Transform a position vector from global inertial reference frame to local
 *  reference frame attached to the turbulence grid.
 */
void global_to_local(const SynthTurbData& turbGrid, const amrex::Real* inp, amrex::Real* out)
{
  const auto* trMat = turbGrid.trMat_;
  amrex::Real in[SynthTurbTraits::NDimMax];
  for (int i=0; i < SynthTurbTraits::NDimMax; ++i)
    in[i] = inp[i] - turbGrid.origin_[i];

  out[0] = trMat[0][0] * in[0] + trMat[0][1] * in[1] + trMat[0][2] * in[2];
  out[1] = trMat[1][0] * in[0] + trMat[1][1] * in[1] + trMat[1][2] * in[2];
  out[2] = trMat[2][0] * in[0] + trMat[2][1] * in[1] + trMat[2][2] * in[2];
}

/** Transform a vector from local reference frame back to the global inertial frame
 *
 */
void local_to_global_vel(const SynthTurbData& turbGrid, const amrex::Real* in, amrex::Real* out)
{
  const auto* trMat = turbGrid.trMat_;
  out[0] = trMat[0][0] * in[0] + trMat[1][0] * in[1] + trMat[2][0] * in[2];
  out[1] = trMat[0][1] * in[0] + trMat[1][1] * in[1] + trMat[2][1] * in[2];
  out[2] = trMat[0][2] * in[0] + trMat[1][2] * in[1] + trMat[2][2] * in[2];
}

/** Determine the left/right indices for a given point along a particular direction
 *
 *  @param turbGrid Turbulence box data
 *  @param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  @param xin Coordinate value in local coordinate frame corresponding to direction provided
 *  @param il Index of the lower bound (populated by this function)
 *  @param ir Index of the upper bound (populated by this function)
 */
void get_lr_indices(
  const SynthTurbData& turbGrid,
  const int dir,
  const amrex::Real xin,
  int& il,
  int& ir)
{
  const amrex::Real xbox =
      xin - std::floor(xin / turbGrid.boxLen_[dir]) * turbGrid.boxLen_[dir];

  il = static_cast<int>(std::floor(xbox / turbGrid.dx_[dir]));
  ir = il + 1;
  if (ir >= turbGrid.boxDims_[dir])
    ir -= turbGrid.boxDims_[dir];
}

/** Determine the left/right indices for a given point along a particular direction
 *
 *  This overload also populates the fractions of left/right states to be used
 *  for interpolations.
 *
 *  @param turbGrid Turbulence box data
 *  @param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  @param xin Coordinate value in local coordinate frame corresponding to
 *  direction provided
 *  @param il Index of the lower bound (populated by this function)
 *  @param ir Index of the upper bound (populated by this function)
 */
void get_lr_indices(
  const SynthTurbData& turbGrid,
  const int dir,
  const amrex::Real xin,
  int& il, int& ir,
  amrex::Real& rxl, amrex::Real& rxr)
{
  const amrex::Real xbox =
      xin - std::floor(xin / turbGrid.boxLen_[dir]) * turbGrid.boxLen_[dir];

  il = static_cast<int>(std::floor(xbox / turbGrid.dx_[dir]));
  ir = il + 1;
  if (ir >= turbGrid.boxDims_[dir])
    ir -= turbGrid.boxDims_[dir];

  const amrex::Real xFrac = xbox - turbGrid.dx_[dir] * il;
  rxl = xFrac / turbGrid.dx_[dir];
  rxr = (1.0 - rxl);
}

/** Indices and interpolation weights for a given point located within the
 *  turbulence box
 */
struct InterpWeights
{
  int il, ir, jl, jr, kl, kr;
  amrex::Real xl, xr, yl, yr, zl, zr;
};

/** Determine if a given point (in local frame) is within the turbulence box
 *
 *  If the point is found within the box, also determine the indices and
 *  interpolation weights for the y and z directions.
 *
 *  @return True if the point is inside the 2-D box
 */
bool find_point_in_box(
  const SynthTurbData& tGrid,
  const amrex::Real* pt,
  InterpWeights& wt)
{
  // Get y and z w.r.t. the lower corner of the grid
  const amrex::Real yy = pt[1] + tGrid.boxLen_[1] * 0.5;
  const amrex::Real zz = pt[2] + tGrid.boxLen_[2] * 0.5;

  bool inBox =
    ((yy >= 0.0) &&
     (yy <= tGrid.boxLen_[1]) &&
     (zz >= 0.0) &&
     (zz <= tGrid.boxLen_[2]));

  if (inBox) {
    get_lr_indices(tGrid, 1, yy, wt.jl, wt.jr, wt.yl, wt.yr);
    get_lr_indices(tGrid, 2, zz, wt.kl, wt.kr, wt.zl, wt.zr);
  }

  return inBox;
}

/** Interpolate the perturbation velocity to a given point from the grid data
 */
void interp_perturb_vel(
  const SynthTurbData& tGrid,
  const InterpWeights& wt,
  amrex::Real* vel)
{
  const int nz = tGrid.boxDims_[2];
  const int nynz = tGrid.boxDims_[1] * tGrid.boxDims_[2];
  // Indices of the 2-D cell that contains the sampling point
  int qidx[4]{wt.jl * nz + wt.kl,
      wt.jr * nz + wt.kl,
      wt.jr * nz + wt.kr,
      wt.jl * nz + wt.kl};

  amrex::Real velL[SynthTurbTraits::NDimMax], velR[SynthTurbTraits::NDimMax];

  // Left quad (t = t)
  velL[0] =
    wt.yl * wt.zl * tGrid.uvel_[qidx[0]] + wt.yr * wt.zl * tGrid.uvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.uvel_[qidx[2]] + wt.yl * wt.zr * tGrid.uvel_[qidx[3]];
  velL[1] =
    wt.yl * wt.zl * tGrid.vvel_[qidx[0]] + wt.yr * wt.zl * tGrid.vvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.vvel_[qidx[2]] + wt.yl * wt.zr * tGrid.vvel_[qidx[3]];
  velL[2] =
    wt.yl * wt.zl * tGrid.wvel_[qidx[0]] + wt.yr * wt.zl * tGrid.wvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.wvel_[qidx[2]] + wt.yl * wt.zr * tGrid.wvel_[qidx[3]];

  for (int i=0; i < 4; ++i)
    qidx[i] += nynz;

  // Right quad (t = t+deltaT)
  velR[0] =
    wt.yl * wt.zl * tGrid.uvel_[qidx[0]] + wt.yr * wt.zl * tGrid.uvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.uvel_[qidx[2]] + wt.yl * wt.zr * tGrid.uvel_[qidx[3]];
  velR[1] =
    wt.yl * wt.zl * tGrid.vvel_[qidx[0]] + wt.yr * wt.zl * tGrid.vvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.vvel_[qidx[2]] + wt.yl * wt.zr * tGrid.vvel_[qidx[3]];
  velR[2] =
    wt.yl * wt.zl * tGrid.wvel_[qidx[0]] + wt.yr * wt.zl * tGrid.wvel_[qidx[1]] +
    wt.yr * wt.zr * tGrid.wvel_[qidx[2]] + wt.yl * wt.zr * tGrid.wvel_[qidx[3]];

  // Interpolation in time
  for (int i=0; i < SynthTurbTraits::NDimMax; ++i)
    vel[i] = wt.xl * velL[i] + wt.xr * velR[i];
}


}

SyntheticTurbulence::SyntheticTurbulence(
  const CFDSim& sim)
    : m_time(sim.time()),
      m_repo(sim.repo()),
      m_mesh(sim.mesh()),
      m_velocity(sim.repo().get_field("velocity")),
      m_density(sim.repo().get_field("density")),
      m_turb_force(sim.repo().declare_field("synth_turb_forcing", 3))
{
  const amrex::Real pi = std::acos(-1.0);

  amrex::ParmParse pp("SynthTurb");

  // NetCDF file containing the turbulence data
  pp.query("turbulence_file", turbFile_.filename);
  process_nc_file(turbFile_, m_turb_grid);

  // Load position and orientation of the grid
  amrex::Real wind_direction;
  pp.query("wind_direction",wind_direction);
  // Convert to radians
  wind_direction *= pi / 180.0;
  amrex::Vector<amrex::Real> location{{0.0,0.0,0.0}};
  pp.queryarr("grid_location",location);

  std::string mean_wind_type = "uniform";
  amrex::Real wind_speed;
  // Default reference height is the center of the turbulence grid
  amrex::Real ref_height = location[2];
  pp.query("mean_wind_speed", wind_speed);
  pp.query("mean_wind_type", mean_wind_type);
  pp.query("mean_wind_ref_height", ref_height);

  if (mean_wind_type == "constant") {
      m_wind_profile.reset(new synth_turb::MeanProfile(wind_speed, ref_height));
  } else if (mean_wind_type == "linear_shear") {
      amrex::Real shear_slope;
      pp.query("shear_slope",shear_slope);
      amrex::Real shear_width;
      pp.query("shear_width",shear_width);
      m_wind_profile.reset(new synth_turb::LinearShearProfile(
                             wind_speed, ref_height, shear_slope, shear_width));
  } else if (mean_wind_type == "power_law") {
      amrex::Real alpha;
      pp.query("power_law_coefficient",alpha);
      m_wind_profile.reset(new synth_turb::PowerLawProfile(wind_speed, ref_height,
                                                         alpha));
  } else {
    throw std::runtime_error("SyntheticTurbulence: invalid mean wind type specified = " + mean_wind_type);
  }

  // Smearing factors
  pp.query("grid_spacing", m_grid_spacing);
  m_epsilon = 2.0 * m_grid_spacing;
  pp.query("gauss_smearing_factor", m_epsilon);
  m_gauss_scaling = 1.0 / (m_epsilon * std::sqrt(pi));

  // Time offsets if any...
  pp.query("time_offset", m_time_offset);

  // Done reading user inputs, process derived data

  // Center of the grid
  m_turb_grid.origin_[0] = location[0];
  m_turb_grid.origin_[1] = location[1];
  m_turb_grid.origin_[2] = location[2];

  // Compute box-fixed reference frame.
  //
  // x-direction points to flow direction (convert from compass direction to vector)
  m_turb_grid.trMat_[0][0] = -std::sin(wind_direction);
  m_turb_grid.trMat_[0][1] = -std::cos(wind_direction);
  m_turb_grid.trMat_[0][2] = 0.0;
  // z always points upwards (for now...)
  m_turb_grid.trMat_[2][0] = 0.0;
  m_turb_grid.trMat_[2][1] = 0.0;
  m_turb_grid.trMat_[2][2] = 1.0;
  // y = z .cross. x
  utils::cross_prod(&m_turb_grid.trMat_[2][0], &m_turb_grid.trMat_[0][0], &m_turb_grid.trMat_[1][0]);

  amrex::Print()
    << "Synthethic turbulence forcing initialized \n"
    << "  Turbulence file = " << turbFile_.filename << "\n"
    << "  Box lengths = [" << m_turb_grid.boxLen_[0] << ", "
    << m_turb_grid.boxLen_[1] << ", " << m_turb_grid.boxLen_[2] << "]\n"
    << "  Box dims = [" << m_turb_grid.boxDims_[0] << ", "
    << m_turb_grid.boxDims_[1] << ", " << m_turb_grid.boxDims_[2] << "]\n"
    << "  Grid dx = [" << m_turb_grid.dx_[0] << ", " << m_turb_grid.dx_[1] << ", "
    << m_turb_grid.dx_[2] << "]\n"
    << "  Centroid (forcing plane) = [" << m_turb_grid.origin_[0] << ", "
    << m_turb_grid.origin_[1] << ", " << m_turb_grid.origin_[2] << "]\n"
    << "  Mean wind profile: U = " << m_wind_profile->reference_velocity()
    << " m/s; Dir = " << wind_direction * 180.0 / pi
    << " deg; H = " << m_wind_profile->reference_height()
    << " m; type = " << mean_wind_type << std::endl;
}

void SyntheticTurbulence::initialize_fields(int level, const amrex::Geometry& geom)
{
    //TODO: Figure out what goes here
}


void SyntheticTurbulence::pre_advance_work()
{
  if (m_is_init)
    initialize();

  update();
}

void SyntheticTurbulence::initialize()
{
  // Convert current time to an equivalent length based on the reference
  // velocity to determine the position within the turbulence grid
  const amrex::Real curTime = m_time.new_time() - m_time_offset;
  const amrex::Real eqivLen = m_wind_profile->reference_velocity() * curTime;
  int il, ir;
  get_lr_indices(m_turb_grid, 0, eqivLen, il, ir);
  load_turb_plane_data(turbFile_, m_turb_grid, il, ir);

  m_is_init = false;
}

void SyntheticTurbulence::update()
{
  // Convert current time to an equivalent length based on the reference
  // velocity to determine the position within the turbulence grid
  const amrex::Real curTime = m_time.new_time() - m_time_offset;
  const amrex::Real eqivLen = m_wind_profile->reference_velocity() * curTime;

  InterpWeights weights;
  get_lr_indices(m_turb_grid, 0, eqivLen,
                 weights.il, weights.ir, weights.xl, weights.xr);

  // Check if we need to refresh the planes
  if (weights.il != m_turb_grid.iLeft_)
    load_turb_plane_data(turbFile_, m_turb_grid, weights.il, weights.ir);


  auto& repo = m_turb_force.repo();
  auto& geom_vec = repo.mesh().Geom();

  const int nlevels = repo.num_active_levels();
  for (int lev=0; lev < nlevels; ++lev) {
      const auto& geom = geom_vec[lev];
      const auto& problo = geom.ProbLoArray();

      const amrex::Real dx = geom.CellSize()[0];
      const amrex::Real dy = geom.CellSize()[1];
      const amrex::Real dz = geom.CellSize()[2];
      const amrex::Real ds = std::cbrt(dx * dy * dz);

      for (amrex::MFIter mfi(m_turb_force(lev)); mfi.isValid(); ++mfi) {
          const auto& bx = mfi.tilebox();
          const auto& turb_force_arr = m_turb_force(lev).array(mfi);
          const auto& vel_arr = m_velocity(lev).array(mfi);
          const auto& rho_arr = m_density(lev).array(mfi);

          amrex::ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            amrex::Real xyzG[SynthTurbTraits::NDimMax];
            amrex::Real xyzL[SynthTurbTraits::NDimMax];
            // velocity in local frame
            amrex::Real velL[SynthTurbTraits::NDimMax];
            // velocity in global frame
            amrex::Real velG[SynthTurbTraits::NDimMax];

            xyzG[0] = problo[0] + (i + 0.5) * dx;
            xyzG[1] = problo[1] + (j + 0.5) * dy;
            xyzG[2] = problo[2] + (k + 0.5) * dz;

            // Transform to local coordinates
            global_to_local(m_turb_grid, xyzG, xyzL);

            InterpWeights wts_loc = weights;

            // Check if the point is in the box, if not we skip this
            // node. The function will also populate the interpolation
            // weights for points that are determined to be within the
            // box.
            bool ptInBox = find_point_in_box(m_turb_grid, xyzL, wts_loc);
            if (ptInBox) {
                // Interpolate perturbation velocities in the local
                // reference frame
                interp_perturb_vel(m_turb_grid, wts_loc, velL);
                // Transform to global coordinates
                local_to_global_vel(m_turb_grid, velL, velG);

                
                // Based on the equations in
                // http://doi.wiley.com/10.1002/we.1608
                // v_n in Eq. 10
                const amrex::Real vMag =
                    std::sqrt(velG[0] * velG[0]
                              + velG[1] * velG[1]
                              + velG[2] * velG[2]);
                // (V_n + 1/2 v_n) in Eq. 10
                const amrex::Real vMagTotal =
                    ((*m_wind_profile)(xyzG[2]) + 0.5 * vMag);
                // Smearing factor (see Eq. 11). The normal direction to
                // the grid is the x-axis of the local reference frame by
                // construction

                const amrex::Real term1 = xyzL[0] / m_epsilon;
                const amrex::Real eta =
                    std::exp(-(term1 * term1)) * m_gauss_scaling;
                const amrex::Real factor =
                    vMagTotal * eta / m_grid_spacing;

                
                turb_force_arr(i,j,k,0) =
                    rho_arr(i,j,k) * velG[0] * factor;
                turb_force_arr(i,j,k,1) =
                    rho_arr(i,j,k) * velG[1] * factor;
                turb_force_arr(i,j,k,2) =
                    rho_arr(i,j,k) * velG[2] * factor;

            }
          });
      }
  }

}

}  // namespace amr_wind
