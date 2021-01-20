#include "amr-wind/physics/SyntheticTurbulence.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

#include "netcdf.h"

namespace amr_wind {
namespace synth_turb {

class LinearShearProfile : public MeanProfile
{
public:
    LinearShearProfile(amrex::Real ref_vel, amrex::Real ref_height,
                       amrex::Real slope, amrex::Real height,
                       int shear_dir)
    : MeanProfile(ref_vel, ref_height, shear_dir),
      m_slope(slope),
      m_half_height(0.5 * height)
  {}

  virtual ~LinearShearProfile() = default;

  virtual amrex::Real operator()(amrex::Real ht) const
  {
    const amrex::Real relHt = ht - m_ref_height;
    if (relHt < -m_half_height)
      return m_ref_vel * (1.0 - m_slope * m_half_height);
    else if (relHt > m_half_height)
      return m_ref_vel * (1.0 + m_slope * m_half_height);
    else
      return m_ref_vel * (1.0 + m_slope * relHt);
  }

private:
  amrex::Real m_slope;
  amrex::Real m_half_height;
};

class PowerLawProfile : public MeanProfile
{
public:
  PowerLawProfile(amrex::Real ref_vel, amrex::Real ref_height, amrex::Real alpha,
                  int shear_dir)
    : MeanProfile(ref_vel, ref_height, shear_dir),
      m_alpha(alpha)
  {}

  virtual ~PowerLawProfile() = default;

  virtual amrex::Real operator()(amrex::Real height) const override
  {
    return m_ref_vel * std::pow((height / m_ref_height), m_alpha);
  }

private:
  const amrex::Real m_alpha;
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
  SyntheticTurbulence::NCBoxTurb& turb_file,
  SynthTurbData& turb_grid)
{
  check_nc_error(nc_open(turb_file.filename.c_str(), NC_NOWRITE, &turb_file.ncid));

  size_t ndim, nx, ny, nz;
  check_nc_error(nc_inq_dimid(turb_file.ncid, "ndim", &turb_file.s_dim));
  check_nc_error(nc_inq_dimlen(turb_file.ncid, turb_file.s_dim, &ndim));
  AMREX_ASSERT(ndim == SynthTurbTraits::n_dim_max);

  // Grid dimensions
  check_nc_error(nc_inq_dimid(turb_file.ncid, "nx", &turb_file.x_dim));
  check_nc_error(nc_inq_dimlen(turb_file.ncid, turb_file.x_dim, &nx));
  check_nc_error(nc_inq_dimid(turb_file.ncid, "ny", &turb_file.y_dim));
  check_nc_error(nc_inq_dimlen(turb_file.ncid, turb_file.y_dim, &ny));
  check_nc_error(nc_inq_dimid(turb_file.ncid, "nz", &turb_file.z_dim));
  check_nc_error(nc_inq_dimlen(turb_file.ncid, turb_file.z_dim, &nz));

  turb_grid.box_dims[0] = nx;
  turb_grid.box_dims[1] = ny;
  turb_grid.box_dims[2] = nz;

  // Box lengths and resolution
  check_nc_error(nc_inq_varid(turb_file.ncid, "box_lengths", &turb_file.boxlen_id));
  check_nc_error(nc_get_var_double(turb_file.ncid, turb_file.boxlen_id, turb_grid.box_len));
  check_nc_error(nc_inq_varid(turb_file.ncid, "dx", &turb_file.dx_id));
  check_nc_error(nc_get_var_double(turb_file.ncid, turb_file.dx_id, turb_grid.dx));

  // Perturbation velocity info
  check_nc_error(nc_inq_varid(turb_file.ncid, "uvel", &turb_file.uid));
  check_nc_error(nc_inq_varid(turb_file.ncid, "vvel", &turb_file.vid));
  check_nc_error(nc_inq_varid(turb_file.ncid, "wvel", &turb_file.wid));
  nc_close(turb_file.ncid);

  // Create data structures to store the perturbation velocities for two planes
  // [t, t+dt] such that the time of interest is within this interval.
  // turb_grid.uvel = SynthTurbTraits::StructField("SynthTurbData::uvel", 2*ny*nz);
  // turb_grid.vvel = SynthTurbTraits::StructField("SynthTurbData::vvel", 2*ny*nz);
  // turb_grid.wvel = SynthTurbTraits::StructField("SynthTurbData::wvel", 2*ny*nz);
  // turb_grid.h_uvel = Kokkos::create_mirror_view(turb_grid.uvel);
  // turb_grid.h_vvel = Kokkos::create_mirror_view(turb_grid.vvel);
  // turb_grid.h_wvel = Kokkos::create_mirror_view(turb_grid.wvel);
  const size_t gridSize = 2 * ny * nz;
  turb_grid.uvel.resize(gridSize);
  turb_grid.vvel.resize(gridSize);
  turb_grid.wvel.resize(gridSize);
}

/** Load two planes of data that bound the current timestep
 *
 *  The data for the y and z directions are loaded for the entire grid at the two planes
 */
void load_turb_plane_data(
  SyntheticTurbulence::NCBoxTurb& turb_file,
  SynthTurbData& turb_grid,
  const int il, const int ir)
{
  check_nc_error(nc_open(turb_file.filename.c_str(), NC_NOWRITE, &turb_file.ncid));

  size_t start[SynthTurbTraits::n_dim_max]{static_cast<size_t>(il), 0, 0};
  size_t count[SynthTurbTraits::n_dim_max]{
    2, static_cast<size_t>(turb_grid.box_dims[1]),
    static_cast<size_t>(turb_grid.box_dims[2])};
  if ((ir - il) == 1) {
    // two consequtive planes load them in one shot
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.uid, start, count, &turb_grid.uvel[0]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.vid, start, count, &turb_grid.vvel[0]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.wid, start, count, &turb_grid.wvel[0]));
  } else {
    // Load the planes separately
    count[0] = 1;
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.uid, start, count, &turb_grid.uvel[0]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.vid, start, count, &turb_grid.vvel[0]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.wid, start, count, &turb_grid.wvel[0]));

    start[0] = static_cast<size_t>(ir);
    const size_t offset = turb_grid.box_dims[1] * turb_grid.box_dims[2];
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.uid, start, count, &turb_grid.uvel[offset]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.vid, start, count, &turb_grid.vvel[offset]));
    check_nc_error(nc_get_vara_double(
      turb_file.ncid, turb_file.wid, start, count, &turb_grid.wvel[offset]));
  }

  // Update left and right indices for future checks
  turb_grid.ileft = il;
  turb_grid.iright = ir;

  nc_close(turb_file.ncid);
}

/** Transform a position vector from global inertial reference frame to local
 *  reference frame attached to the turbulence grid.
 */
void global_to_local(const SynthTurbData& turb_grid, const amrex::Real* inp, amrex::Real* out)
{
  const auto* tr_mat = turb_grid.tr_mat;
  amrex::Real in[SynthTurbTraits::n_dim_max];
  for (int i=0; i < SynthTurbTraits::n_dim_max; ++i)
    in[i] = inp[i] - turb_grid.origin[i];

  out[0] = tr_mat[0][0] * in[0] + tr_mat[0][1] * in[1] + tr_mat[0][2] * in[2];
  out[1] = tr_mat[1][0] * in[0] + tr_mat[1][1] * in[1] + tr_mat[1][2] * in[2];
  out[2] = tr_mat[2][0] * in[0] + tr_mat[2][1] * in[1] + tr_mat[2][2] * in[2];
}

/** Transform a vector from local reference frame back to the global inertial frame
 *
 */
void local_to_global_vel(const SynthTurbData& turb_grid, const amrex::Real* in, amrex::Real* out)
{
  const auto* tr_mat = turb_grid.tr_mat;
  out[0] = tr_mat[0][0] * in[0] + tr_mat[1][0] * in[1] + tr_mat[2][0] * in[2];
  out[1] = tr_mat[0][1] * in[0] + tr_mat[1][1] * in[1] + tr_mat[2][1] * in[2];
  out[2] = tr_mat[0][2] * in[0] + tr_mat[1][2] * in[1] + tr_mat[2][2] * in[2];
}

/** Determine the left/right indices for a given point along a particular direction
 *
 *  @param turb_grid Turbulence box data
 *  @param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  @param xin Coordinate value in local coordinate frame corresponding to direction provided
 *  @param il Index of the lower bound (populated by this function)
 *  @param ir Index of the upper bound (populated by this function)
 */
void get_lr_indices(
  const SynthTurbData& turb_grid,
  const int dir,
  const amrex::Real xin,
  int& il,
  int& ir)
{
  const amrex::Real xbox =
      xin - std::floor(xin / turb_grid.box_len[dir]) * turb_grid.box_len[dir];

  il = static_cast<int>(std::floor(xbox / turb_grid.dx[dir]));
  ir = il + 1;
  if (ir >= turb_grid.box_dims[dir])
    ir -= turb_grid.box_dims[dir];
}

/** Determine the left/right indices for a given point along a particular direction
 *
 *  This overload also populates the fractions of left/right states to be used
 *  for interpolations.
 *
 *  @param turb_grid Turbulence box data
 *  @param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  @param xin Coordinate value in local coordinate frame corresponding to
 *  direction provided
 *  @param il Index of the lower bound (populated by this function)
 *  @param ir Index of the upper bound (populated by this function)
 */
void get_lr_indices(
  const SynthTurbData& turb_grid,
  const int dir,
  const amrex::Real xin,
  int& il, int& ir,
  amrex::Real& rxl, amrex::Real& rxr)
{
  const amrex::Real xbox =
      xin - std::floor(xin / turb_grid.box_len[dir]) * turb_grid.box_len[dir];

  il = static_cast<int>(std::floor(xbox / turb_grid.dx[dir]));
  ir = il + 1;
  if (ir >= turb_grid.box_dims[dir])
    ir -= turb_grid.box_dims[dir];

  const amrex::Real xFrac = xbox - turb_grid.dx[dir] * il;
  rxl = xFrac / turb_grid.dx[dir];
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
  const SynthTurbData& t_grid,
  const amrex::Real* pt,
  InterpWeights& wt)
{
  // Get y and z w.r.t. the lower corner of the grid
  const amrex::Real yy = pt[1] + t_grid.box_len[1] * 0.5;
  const amrex::Real zz = pt[2] + t_grid.box_len[2] * 0.5;

  bool inBox =
    ((yy >= 0.0) &&
     (yy <= t_grid.box_len[1]) &&
     (zz >= 0.0) &&
     (zz <= t_grid.box_len[2]));

  if (inBox) {
    get_lr_indices(t_grid, 1, yy, wt.jl, wt.jr, wt.yl, wt.yr);
    get_lr_indices(t_grid, 2, zz, wt.kl, wt.kr, wt.zl, wt.zr);
  }

  return inBox;
}

/** Interpolate the perturbation velocity to a given point from the grid data
 */
void interp_perturb_vel(
  const SynthTurbData& t_grid,
  const InterpWeights& wt,
  amrex::Real* vel)
{
  const int nz = t_grid.box_dims[2];
  const int nynz = t_grid.box_dims[1] * t_grid.box_dims[2];
  // Indices of the 2-D cell that contains the sampling point
  int qidx[4]{wt.jl * nz + wt.kl,
      wt.jr * nz + wt.kl,
      wt.jr * nz + wt.kr,
      wt.jl * nz + wt.kl};

  amrex::Real vel_l[SynthTurbTraits::n_dim_max], vel_r[SynthTurbTraits::n_dim_max];

  // Left quad (t = t)
  vel_l[0] =
    wt.yl * wt.zl * t_grid.uvel[qidx[0]] + wt.yr * wt.zl * t_grid.uvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.uvel[qidx[2]] + wt.yl * wt.zr * t_grid.uvel[qidx[3]];
  vel_l[1] =
    wt.yl * wt.zl * t_grid.vvel[qidx[0]] + wt.yr * wt.zl * t_grid.vvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.vvel[qidx[2]] + wt.yl * wt.zr * t_grid.vvel[qidx[3]];
  vel_l[2] =
    wt.yl * wt.zl * t_grid.wvel[qidx[0]] + wt.yr * wt.zl * t_grid.wvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.wvel[qidx[2]] + wt.yl * wt.zr * t_grid.wvel[qidx[3]];

  for (int i=0; i < 4; ++i)
    qidx[i] += nynz;

  // Right quad (t = t+deltaT)
  vel_r[0] =
    wt.yl * wt.zl * t_grid.uvel[qidx[0]] + wt.yr * wt.zl * t_grid.uvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.uvel[qidx[2]] + wt.yl * wt.zr * t_grid.uvel[qidx[3]];
  vel_r[1] =
    wt.yl * wt.zl * t_grid.vvel[qidx[0]] + wt.yr * wt.zl * t_grid.vvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.vvel[qidx[2]] + wt.yl * wt.zr * t_grid.vvel[qidx[3]];
  vel_r[2] =
    wt.yl * wt.zl * t_grid.wvel[qidx[0]] + wt.yr * wt.zl * t_grid.wvel[qidx[1]] +
    wt.yr * wt.zr * t_grid.wvel[qidx[2]] + wt.yl * wt.zr * t_grid.wvel[qidx[3]];

  // Interpolation in time
  for (int i=0; i < SynthTurbTraits::n_dim_max; ++i)
    vel[i] = wt.xl * vel_l[i] + wt.xr * vel_r[i];
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
  pp.query("turbulence_file", m_turb_file.filename);
  process_nc_file(m_turb_file, m_turb_grid);

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
      int shear_dir=2;
      pp.query("shear_dir",shear_dir);
      m_wind_profile.reset(new synth_turb::LinearShearProfile(
                               wind_speed, ref_height, shear_slope,
                               shear_width, shear_dir));
  } else if (mean_wind_type == "power_law") {
      amrex::Real alpha;
      pp.query("power_law_coefficient",alpha);
      int shear_dir=2;
      pp.query("shear_dir",shear_dir);
      m_wind_profile.reset(new synth_turb::PowerLawProfile(wind_speed,
                                                           ref_height,
                                                           alpha, shear_dir));
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
  m_turb_grid.origin[0] = location[0];
  m_turb_grid.origin[1] = location[1];
  m_turb_grid.origin[2] = location[2];

  // Compute box-fixed reference frame.
  //
  // x-direction points to flow direction (convert from compass direction to vector)
  m_turb_grid.tr_mat[0][0] = -std::sin(wind_direction);
  m_turb_grid.tr_mat[0][1] = -std::cos(wind_direction);
  m_turb_grid.tr_mat[0][2] = 0.0;
  // z always points upwards (for now...)
  m_turb_grid.tr_mat[2][0] = 0.0;
  m_turb_grid.tr_mat[2][1] = 0.0;
  m_turb_grid.tr_mat[2][2] = 1.0;
  // y = z .cross. x
  utils::cross_prod(&m_turb_grid.tr_mat[2][0], &m_turb_grid.tr_mat[0][0], &m_turb_grid.tr_mat[1][0]);

  amrex::Print()
    << "Synthethic turbulence forcing initialized \n"
    << "  Turbulence file = " << m_turb_file.filename << "\n"
    << "  Box lengths = [" << m_turb_grid.box_len[0] << ", "
    << m_turb_grid.box_len[1] << ", " << m_turb_grid.box_len[2] << "]\n"
    << "  Box dims = [" << m_turb_grid.box_dims[0] << ", "
    << m_turb_grid.box_dims[1] << ", " << m_turb_grid.box_dims[2] << "]\n"
    << "  Grid dx = [" << m_turb_grid.dx[0] << ", " << m_turb_grid.dx[1] << ", "
    << m_turb_grid.dx[2] << "]\n"
    << "  Centroid (forcing plane) = [" << m_turb_grid.origin[0] << ", "
    << m_turb_grid.origin[1] << ", " << m_turb_grid.origin[2] << "]\n"
    << "  Mean wind profile: U = " << m_wind_profile->reference_velocity()
    << " m/s; Dir = " << wind_direction * 180.0 / pi
    << " deg; H = " << m_wind_profile->reference_height()
    << " m; type = " << mean_wind_type << std::endl;
}

void SyntheticTurbulence::initialize_fields(
    int // level
    ,
    const amrex::Geometry& // geom)
)
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
  load_turb_plane_data(m_turb_file, m_turb_grid, il, ir);

  m_is_init = false;
}

void SyntheticTurbulence::update()
{
  // Convert current time to an equivalent length based on the reference
  // velocity to determine the position within the turbulence grid
  const amrex::Real cur_time = m_time.new_time() - m_time_offset;
  const amrex::Real eqiv_len = m_wind_profile->reference_velocity() * cur_time;

  InterpWeights weights;
  get_lr_indices(m_turb_grid, 0, eqiv_len,
                 weights.il, weights.ir, weights.xl, weights.xr);

  // Check if we need to refresh the planes
  if (weights.il != m_turb_grid.ileft)
    load_turb_plane_data(m_turb_file, m_turb_grid, weights.il, weights.ir);

  auto& repo = m_turb_force.repo();
  auto& geom_vec = repo.mesh().Geom();

  const int sdir = (*m_wind_profile).shear_dir();

  const int nlevels = repo.num_active_levels();
  for (int lev=0; lev < nlevels; ++lev) {
      const auto& geom = geom_vec[lev];
      const auto& problo = geom.ProbLoArray();

      const amrex::Real dx = geom.CellSize()[0];
      const amrex::Real dy = geom.CellSize()[1];
      const amrex::Real dz = geom.CellSize()[2];

      for (amrex::MFIter mfi(m_turb_force(lev)); mfi.isValid(); ++mfi) {
          const auto& bx = mfi.tilebox();
          const auto& turb_force_arr = m_turb_force(lev).array(mfi);
          const auto& rho_arr = m_density(lev).array(mfi);

          amrex::ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            amrex::Real xyz_g[SynthTurbTraits::n_dim_max];
            amrex::Real xyz_l[SynthTurbTraits::n_dim_max];
            // velocity in local frame
            amrex::Real vel_l[SynthTurbTraits::n_dim_max];
            // velocity in global frame
            amrex::Real vel_g[SynthTurbTraits::n_dim_max];

            xyz_g[0] = problo[0] + (i + 0.5) * dx;
            xyz_g[1] = problo[1] + (j + 0.5) * dy;
            xyz_g[2] = problo[2] + (k + 0.5) * dz;

            // Transform to local coordinates
            global_to_local(m_turb_grid, xyz_g, xyz_l);

            InterpWeights wts_loc = weights;

            // Check if the point is in the box, if not we skip this
            // node. The function will also populate the interpolation
            // weights for points that are determined to be within the
            // box.
            bool ptInBox = find_point_in_box(m_turb_grid, xyz_l, wts_loc);
            if (ptInBox) {
                // Interpolate perturbation velocities in the local
                // reference frame
                interp_perturb_vel(m_turb_grid, wts_loc, vel_l);
                // Transform to global coordinates
                local_to_global_vel(m_turb_grid, vel_l, vel_g);

                // Based on the equations in
                // http://doi.wiley.com/10.1002/we.1608
                // v_n in Eq. 10
                const amrex::Real v_mag =
                    std::sqrt(vel_g[0] * vel_g[0]
                              + vel_g[1] * vel_g[1]
                              + vel_g[2] * vel_g[2]);
                // (V_n + 1/2 v_n) in Eq. 10
                const amrex::Real v_mag_total =
                    ((*m_wind_profile)(xyz_g[sdir]) + 0.5 * v_mag);
                // Smearing factor (see Eq. 11). The normal direction to
                // the grid is the x-axis of the local reference frame by
                // construction

                const amrex::Real term1 = xyz_l[0] / m_epsilon;
                const amrex::Real eta =
                    std::exp(-(term1 * term1)) * m_gauss_scaling;
                const amrex::Real factor =
                    v_mag_total * eta;

                turb_force_arr(i,j,k,0) =
                    rho_arr(i,j,k) * vel_g[0] * factor;
                turb_force_arr(i,j,k,1) =
                    rho_arr(i,j,k) * vel_g[1] * factor;
                turb_force_arr(i,j,k,2) =
                    rho_arr(i,j,k) * vel_g[2] * factor;

            }
          });
      }
  }

}

}  // namespace amr_wind
