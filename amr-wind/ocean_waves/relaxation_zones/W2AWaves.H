#ifndef W2AWAVES_H
#define W2AWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#include <vector>
#include <complex>
#ifdef AMR_WIND_USE_W2A
#include "read_modes.h"
#include "modes_hosgrid.h"
#include "interp_to_mfab.h"
#endif

namespace amr_wind::ocean_waves {

struct W2AWavesData : public RelaxZonesBaseData
{
    // Prefix for HOS files
    std::string modes_file{"modes_HOS_SWENSE.dat"};
    // Number of modes in each direction
    int n0{0};
    int n1{0};
    // Spatial resolution in each direction
    amrex::Real dx0{0.};
    amrex::Real dx1{0.};
    // Depth of reference sim
    amrex::Real depth{0.};
    // Length to dimensionalize nondim data
    amrex::Real dimL{0.};
    // Timestep of reference data
    amrex::Real dt_modes{0.0};
    // Time when ow arrays were last changed (last timestep)
    amrex::Real t_last{-1.0};
    // Time of current HOS lists
    amrex::Real t{0.0};
    // Timestep of current HOS lists
    int ntime{0};

    // Vectors to store modes
    std::vector<std::complex<double>> mX, mY, mZ, mFS;

// Struct variables that have special types
#ifdef AMR_WIND_USE_W2A
    // FFTW plan
    fftw_plan plan;
    // FFTW pointers to store modes ready for ifft
    fftw_complex *eta_mptr, *u_mptr, *v_mptr, *w_mptr;
    // ReadModes object
    ReadModes rmodes;
#endif

    // Height vector (where velocity is sampled) stuff
    amrex::Vector<amrex::Real> hvec;
    // Cell height of interpolation grid near interface
    amrex::Real dz0{0.};
    // Number of sampling points
    int nheights{0};
    // Vectors for spatial data from transformed modes
    amrex::Gpu::DeviceVector<amrex::Real> sp_eta_vec, sp_u_vec, sp_v_vec,
        sp_w_vec;
    // Spacing of spatial data
    amrex::Real sp_dx{0.}, sp_dy{0.};
};

struct W2AWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = W2AWavesData;
    using DataType = OceanWavesDataHolder<W2AWaves>;

    static std::string identifier() { return "W2AWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* W2AWAVES_H */
