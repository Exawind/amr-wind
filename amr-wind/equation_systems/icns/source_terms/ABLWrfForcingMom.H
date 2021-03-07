#ifndef ABLWRFFORCINGMOM
#define ABLWRFFORCINGMOM

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Forcing term for meso-micro scale coupled ABL simulations
 *  \ingroup icns_src we_abl
 *
 *  \sa ABL
 */

class WRFfile{

 public:

  explicit WRFfile();

  ~WRFfile();

  void read_forcing_file(){

    auto ncf = ncutils::NCFile::open_par(
        m_wrf_file, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    m_nheight  = ncf.dim("nheight").len();
    m_ntime = ncf.dim("ntime").len();

    m_wrf_height.resize(m_nheight);
    m_wrf_time.resize(m_ntime);

    ncf.var("heights").get(m_wrf_height.data());
    ncf.var("times").get(m_wrf_time.data());

    m_wrf_mom.resize(m_nheight*m_ntime*2);
    m_wrf_temp.resize(m_nheight*m_ntime);
    m_wrf_tflux.resize(m_ntime); 
    
    ncf.var("wrf_momentum").get(m_wrf_mom.data());
    ncf.var("wrf_temperature").get(m_wrf_temp.data());
    ncf.var("wrf_tflux").get(m_wrf_tflux.data());

  };

  void set_forcing_file(std::string filename) {
    m_wrf_file = filename;
  };
  
 private: 
  std::string m_wrf_file;

  amrex::Vector<amrex::Real> m_wrf_height;
  amrex::Vector<amrex::Real> m_wrf_time;

  amrex::Vector<amrex::Real> m_wrf_mom;
  amrex::Vector<amrex::Real> m_wrf_temp;
  amrex::Vector<amrex::Real> m_wrf_tflux;

  int m_nheight;
  int m_ntime;

};

class ABLWrfForcingMom: public MomentumSource::Register<ABLWrfForcingMom>, public WRFfile
{
public:

  static const std::string identifier() { return "ABLWrfForcingMom"; }
  
  explicit ABLWrfForcingMom(const CFDSim& sim);
  
  virtual ~ABLWrfForcingMom();

  virtual void operator()(
      const int lev,
      const amrex::MFIter& mfi,
      const amrex::Box& bx,
      const FieldState fstate,
      const amrex::Array4<amrex::Real>& src_term) const override;

  void mean_velocity_heights(VelPlaneAveraging&);

  void mean_velocity_init(VelPlaneAveraging&);
 

 private:
  const SimTime& m_time;

  const amrex::AmrCore& m_mesh;
  
  std::string m_forcing_scheme{"direct"};

  amrex::Gpu::DeviceVector<amrex::Real> m_velAvg_ht;
  amrex::Gpu::DeviceVector<amrex::Real> m_uAvg_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_vAvg_vals;

  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_u_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_v_vals;
  
  int m_idx_time;

  int m_axis{2};

};

} // namespace icns
} // namespace pde
} // namespace amr_wind


#endif // ABLWRFFORCINGMOM_H

