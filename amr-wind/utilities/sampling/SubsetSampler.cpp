#include "amr-wind/utilities/sampling/SubsetSampler.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

#include "amr-wind/core/Field.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace sampling {

SubsetSampler::SubsetSampler(const CFDSim& sim)
  : m_sim(sim)
  , m_mesh(sim.mesh())
  , m_velocity(sim.repo().get_field("velocity"))
  , m_mesh_fac_cc(sim.repo().get_field("mesh_scaling_factor_cc"))
  , m_nu_coord_cc(sim.repo().get_field("non_uniform_coord_cc"))
{}

SubsetSampler::~SubsetSampler() = default;

void SubsetSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    
    // Get the input parameters
    pp.getarr("xlim", m_xlim);
    pp.getarr("ylim", m_ylim);
    pp.getarr("zlim", m_zlim);
    pp.query("level", m_level);

    // Check the inputs
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_xlim.size()) == 2);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_ylim.size()) == 2);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_zlim.size()) == 2);
    
    AMREX_ALWAYS_ASSERT(m_xlim[1] >= m_xlim[0]);
    AMREX_ALWAYS_ASSERT(m_ylim[1] >= m_ylim[0]);
    AMREX_ALWAYS_ASSERT(m_zlim[1] >= m_zlim[0]);

    const int level   = m_level;

    auto& velocity    = m_velocity(level);
    auto& nu_coord_cc = m_nu_coord_cc(level);

    const int Nproc   = amrex::ParallelDescriptor::NProcs();
    const int iproc   = amrex::ParallelDescriptor::MyProc();

    amrex::Vector<int> proc_cc_count(Nproc, 0);

    int npt=0;
    // Loop through and count total number of points inside box
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
	const auto& nu_cc = nu_coord_cc.array(mfi);
	amrex::Loop(vbx, [=, &npt](int i, int j, int k) noexcept {
	                amrex::Real x = nu_cc(i, j, k, 0);
			amrex::Real y = nu_cc(i, j, k, 1);
			amrex::Real z = nu_cc(i, j, k, 2);
			if (   (m_xlim[0] <= x) && (x <= m_xlim[1]) 
			    && (m_ylim[0] <= y) && (y <= m_ylim[1])
			    && (m_zlim[0] <= z) && (z <= m_zlim[1])) {
			  npt++;
			}
	  });
    }

    // Gather and sum the total point counts
#ifdef AMREX_USE_MPI
    MPI_Allgather(&npt, 1, MPI_INT, proc_cc_count.data(), 1, MPI_INT,
		  amrex::ParallelDescriptor::Communicator());
#else
    proc_cc_count[0] = npt;
#endif
    m_npts = std::accumulate(proc_cc_count.begin(), proc_cc_count.end(), 0);

    // Compute the offsets for storage
    m_proc_offsets.resize(Nproc);
    m_proc_offsets[0] = 0;
    for (int i=1; i<Nproc; i++)
      m_proc_offsets[i] = m_proc_offsets[i-1] + proc_cc_count[i-1];

    // Go through and get the coordinates of the points inside the box
    const auto& dx      = m_mesh.Geom(level).CellSizeArray();
    const auto* prob_lo = m_sim.mesh().Geom(level).ProbLo();

    m_pos_nu_vec.resize(m_npts, vs::Vector::zero());
    m_pos_comp_vec.resize(m_npts, vs::Vector::zero());
							 
    auto *pos_nu_arr   = m_pos_nu_vec.data();
    auto *pos_comp_arr = m_pos_comp_vec.data();

    int idx = 0;
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
	const auto& nu_cc = nu_coord_cc.array(mfi);
	const int offset  = m_proc_offsets[iproc];
	amrex::Loop(vbx, [=, &idx] (int i, int j, int k) 
		    noexcept {
	                amrex::Real x = nu_cc(i, j, k, 0);
			amrex::Real y = nu_cc(i, j, k, 1);
			amrex::Real z = nu_cc(i, j, k, 2);
			if (   (m_xlim[0] <= x) && (x <= m_xlim[1]) 
			    && (m_ylim[0] <= y) && (y <= m_ylim[1])
			    && (m_zlim[0] <= z) && (z <= m_zlim[1])) {
			  // Add points to local array
			  auto &pos_nu   = pos_nu_arr[idx+offset];
			  auto &pos_comp = pos_comp_arr[idx+offset];
			  for (int n=0; n<AMREX_SPACEDIM; n++){
			    pos_nu[n] = nu_cc(i, j, k, n);
			  }
			  pos_comp[0] = prob_lo[0] + (i + 0.5)*dx[0];
			  pos_comp[1] = prob_lo[1] + (j + 0.5)*dx[1];
			  pos_comp[2] = prob_lo[2] + (k + 0.5)*dx[2];

			  // advance to next point
			  idx++;
			}
	  });
    }

#ifdef AMREX_USE_MPI
    MPI_Allreduce(
        MPI_IN_PLACE, &(m_pos_nu_vec[0][0]), m_npts*AMREX_SPACEDIM,
        MPI_DOUBLE, MPI_SUM, amrex::ParallelDescriptor::Communicator());
    MPI_Allreduce(
        MPI_IN_PLACE, &(m_pos_comp_vec[0][0]), m_npts*AMREX_SPACEDIM,
        MPI_DOUBLE, MPI_SUM, amrex::ParallelDescriptor::Communicator());
#endif

}

void SubsetSampler::sampling_locations(SampleLocType& locs) const
{
    locs.resize(m_npts);
    for (int i=0; i<m_npts; i++) {
      for (int d=0; d < AMREX_SPACEDIM; d++) {
	locs[i][d] = m_pos_comp_vec[i][d];
      }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void SubsetSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("xlim",  m_xlim);
    grp.put_attr("ylim",  m_ylim);
    grp.put_attr("zlim",  m_zlim);
    grp.def_scalar("level", NC_INT);
    grp.def_var("physical_coordinates", NC_DOUBLE, {"num_points", "ndim"});

}
void SubsetSampler::populate_netcdf_metadata(const ncutils::NCGroup& grp) const 
{
    auto grplevel = grp.var("level");
    grplevel.put(&m_level);
  
    auto xyz = grp.var("physical_coordinates");
    std::vector<size_t> start{0, 0};
    std::vector<size_t> count{(size_t)m_npts, AMREX_SPACEDIM};
    xyz.put(&m_pos_nu_vec[0][0], start, count);

}

#else
void SubsetSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void SubsetSampler::populate_netcdf_metadata(const ncutils::NCGroup&) const 
{}
#endif

} // namespace sampling
} // namespace amr_wind
