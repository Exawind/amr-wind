#include "amr-wind/wind_energy/ABLiStats.H"
#include "amr-wind/derive/derive_K.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/DirectionSelector.H"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"

namespace amr_wind {

struct temp_grad{
  double grad_z;
  double max_grad_loc;
};

ABLiStats::ABLiStats(CFDSim& sim, const ABLWallFunction& abl_wall_func)
    : m_sim(sim)
    , m_abl_wall_func(abl_wall_func)
    , m_temperature(sim.repo().get_field("temperature"))
{}

ABLiStats::~ABLiStats() = default;

void ABLiStats::initialize()
{
    BL_PROFILE("amr-wind::ABLiStats::initialize");

    {
        amrex::ParmParse pp("ABL");
        pp.query("stats_output_frequency", m_out_freq);
        pp.query("stats_output_format", m_out_fmt);
        pp.query("normal_direction", m_normal_dir);
        AMREX_ASSERT((0 <= m_normal_dir) && (m_normal_dir < AMREX_SPACEDIM));
    }

    //Get normal direction and associated stuff
    const auto& geom = (this->m_sim.repo()).mesh().Geom()[0];
    amrex::Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    switch (m_normal_dir) {
    case 0:
        m_ncells_h1 = dhi.y-dlo.y+1;
        m_ncells_h2 = dhi.z-dlo.z+1;
        break;
    case 1:
        m_ncells_h1 = dhi.x-dlo.x+1;
        m_ncells_h2 = dhi.z-dlo.z+1;
        break;
    case 2:
        m_ncells_h1 = dhi.x-dlo.x+1;
        m_ncells_h2 = dhi.y-dlo.y+1;
        break;
    }
    m_dn = geom.CellSize()[m_normal_dir];
    
    if (m_out_fmt == "netcdf") prepare_netcdf_file();
}

void ABLiStats::post_advance_work()
{
    BL_PROFILE("amr-wind::ABLiStats::post_advance_work");
    const auto& time = m_sim.time();
    const int tidx = time.time_index();
    // Skip processing if it is not an output timestep
    if (!(tidx % m_out_freq == 0)) return;

    switch (m_normal_dir) {
    case 0:
        compute_zi(YDir(), ZDir());
        break;
    case 1:
        compute_zi(XDir(), ZDir());
        break;
    case 2:
        compute_zi(XDir(), YDir());
        break;
    }

    process_output();
}

template <typename h1_dir, typename h2_dir> 
void ABLiStats::compute_zi(const h1_dir& h1Sel, const h2_dir& h2Sel) {

    auto gradT = (this->m_sim.repo()).create_scratch_field(3,m_temperature.num_grow()[0]);
    compute_gradient(*gradT, m_temperature);

    //Only compute zi using coarsest level
    
    amrex::Gpu::ManagedVector<temp_grad> tgrad(m_ncells_h1 * m_ncells_h2);
    auto * tgrad_ptr = tgrad.data();
    amrex::Gpu::ManagedVector<temp_grad> tgrad_g(m_ncells_h1 * m_ncells_h2);
    auto * tgrad_g_ptr = tgrad_g.data();
    
    for (amrex::MFIter mfi(m_temperature(0)); mfi.isValid(); ++mfi) {
      const auto& bx = mfi.tilebox();
      const auto& gradT_arr = (*gradT)(0).array(mfi);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int h1 = h1Sel(i,j,k);
        int h2 = h2Sel(i,j,k);        
        if (tgrad[h2*m_ncells_h2 + h1].grad_z
            < gradT_arr(i,j,k,m_normal_dir)) {
          tgrad_ptr[h2*m_ncells_h2 + h1].grad_z
              = gradT_arr(i,j,k,m_normal_dir);
          tgrad_ptr[h2*m_ncells_h2 + h1].max_grad_loc = (k+0.5)*m_dn;
        }
      });
    }

    MPI_Reduce(tgrad_ptr, tgrad_g_ptr, m_ncells_h1 * m_ncells_h2,
               MPI_2DOUBLE_PRECISION, MPI_MAXLOC,
               amrex::ParallelDescriptor::IOProcessorNumber(),
               amrex::ParallelDescriptor::Communicator() );

    m_zi = 0.0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
      for (int i=0; i < m_ncells_h1 * m_ncells_h2; i++) {
          m_zi += tgrad_g_ptr[i].max_grad_loc;
      }
      m_zi /= (m_ncells_h1 * m_ncells_h2);
    }
}

void ABLiStats::process_output()
{
    if (m_out_fmt == "ascii") {
        write_ascii();
    } else if (m_out_fmt == "netcdf") {
        write_netcdf();
    } else {
        amrex::Abort("ABLiStats: Invalid output format encountered");
    }
}

void ABLiStats::write_ascii()
{
    BL_PROFILE("amr-wind::ABLiStats::write_ascii");
    amrex::Print() << "WARNING: ABLiStats: ASCII output will impact performance"
                   << std::endl;

    const std::string istat_dir = "i_stats";
    const std::string sname = amrex::Concatenate("i_stats", m_sim.time().time_index());

    if (!amrex::UtilCreateDirectory(istat_dir, 0755)) {
        amrex::CreateDirectoryFailed(istat_dir);
    }
    const std::string fname = istat_dir + "/" + sname + ".txt";
    //TODO: Do something to write the data to file here
}

void ABLiStats::prepare_netcdf_file()
{
#ifdef AMR_WIND_USE_NETCDF

    const std::string istat_dir = "i_stats";
    const std::string sname = amrex::Concatenate("i_stats", m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(istat_dir, 0755)) {
        amrex::CreateDirectoryFailed(istat_dir);
    }
    m_ncfile_name = istat_dir + "/" + sname + ".nc";

    // Only I/O processor handles NetCDF generation
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    auto ncf = ncutils::NCFile::create(m_ncfile_name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string npart_name = "num_points";
    const std::vector<std::string> two_dim{nt_name, npart_name};
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind data sampling output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_DOUBLE, {nt_name});
    ncf.def_var("ustar", NC_DOUBLE, {nt_name});
    ncf.def_var("wstar", NC_DOUBLE, {nt_name});
    ncf.def_var("L", NC_DOUBLE, {nt_name});
    ncf.def_var("zi", NC_DOUBLE, {nt_name});
    ncf.exit_def_mode();

#else
    amrex::Abort(
        "NetCDF support was not enabled during build time. Please recompile or "
        "use native format");
#endif
}

void ABLiStats::write_netcdf()
{
#ifdef AMR_WIND_USE_NETCDF

    if (!amrex::ParallelDescriptor::IOProcessor()) return;
    auto ncf = ncutils::NCFile::open(m_ncfile_name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of the next timestep
    const size_t nt = ncf.dim(nt_name).len();
    {
        auto time = m_sim.time().new_time();
        ncf.var("time").put(&time, {nt}, {1});
        auto ustar = m_abl_wall_func.utau();
        ncf.var("ustar").put(&ustar, {nt}, {1});
        double wstar = 0.0;
        ncf.var("wstar").put(&wstar, {nt}, {1});
        double L = 1e18;
        ncf.var("L").put(&L, {nt}, {1});
        ncf.var("zi").put(&m_zi, {nt}, {1});
    }
    ncf.close();
#endif
}

}
