#ifndef ABLWRFFORCINGMOM
#define ABLWRFFORCINGMOM

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

class PlaneAveraging;

namespace amr_wind {
namespace pde {
namespace icns {

/** Forcing term for meso-micro scale coupled ABL simulations
 *  \ingroup icns_src we_abl
 *
 *  \sa ABL
 */

class ABLWrfForcingMom: public MomentumSource::Register<ABLWrfForcingMom>
{
public:
  explicit ABLWrfForcingMom(){const CFDSim& sim};
  
  virtual ~ABLWrfForcingMom();

  virtual void operator()(
      const int lev,
      const amrex::MFIter& mfi,
      const amrex::Box& bx,
      const FieldState fstate,
      const amrex::Array4<amrex::Real>& src_term) const override;

  void read_forcing_file();

  void mean_velocity_heights(const FieldPlaneAveraging&);
 

 private:
  const SimTime& m_time;

  std::string m_wrf_file;
  std::string m_forcing_scheme{"direct"};

  amrex::Array<amrex::Real> m_wrf_mom;

  amrex::Vector<amrex::Real> m_wrf_height;
  amrex::Vector<amrex::Real> m_wrf_time;

  int m_nheight;
  int m_ntime;

  amrex::Gpu::DeviceVector<amrex::Real> m_velAvg_ht;
  amrex::Gpu::DeviceVector<amrex::Real> m_uAvg_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_vAvg_vals;

  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_u_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_v_vals;

  
  int m_idx_time;

};

} // namespace icns
} // namespace pde
} // namespace amr_wind


#endif // ABLWRFFORCINGMOM_H

