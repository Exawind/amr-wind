#ifndef ABLWRFFORCINGMOM
#define ABLWRFFORCINGMOM

#include <memory>
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/wind_energy/ABLWrf.H"

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

  static const std::string identifier() { return "ABLWrfForcingMom"; }
  
  explicit ABLWrfForcingMom(const CFDSim& sim);
  
  virtual ~ABLWrfForcingMom();

  virtual void operator()(
      const int lev,
      const amrex::MFIter& mfi,
      const amrex::Box& bx,
      const FieldState fstate,
      const amrex::Array4<amrex::Real>& src_term) const override;

  void mean_velocity_heights(const VelPlaneAveraging&, std::unique_ptr<ABLWRFfile>&);

  void mean_velocity_init(const VelPlaneAveraging&, const ABLWRFfile&);

  void invertMat(const amrex::Array2D<amrex::Real,0,3,0,3>&, amrex::Array2D<amrex::Real,0,3,0,3>&);

  void indirectForcingInit();
 
 private:
  const SimTime& m_time;

  const amrex::AmrCore& m_mesh;
  
  std::string m_forcing_scheme;

  amrex::Gpu::DeviceVector<amrex::Real> m_velAvg_ht;
  amrex::Gpu::DeviceVector<amrex::Real> m_uAvg_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_vAvg_vals;

  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_u_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_v_vals;

  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_ht;

  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_avg_error;

  amrex::Gpu::DeviceVector<amrex::Real> m_error_wrf_avg_U;
  amrex::Gpu::DeviceVector<amrex::Real> m_error_wrf_avg_V;

  amrex::Vector<amrex::Real> m_zht;
  
  int m_idx_time;

  int m_axis{2};

  int m_ind_polyOrder;

  int m_nht;

  amrex::Array2D<amrex::Real,0,3,0,3> m_im_zTz;

  amrex::Array<amrex::Real, 4> m_poly_coeff_U;
  amrex::Array<amrex::Real, 4> m_poly_coeff_V;

};

} // namespace icns
} // namespace pde
} // namespace amr_wind


#endif // ABLWRFFORCINGMOM_H

