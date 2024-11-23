#ifndef W2AWAVES_H
#define W2AWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#include <vector>
#include <complex>
#ifdef AMR_WIND_USE_W2A
#include "Waves2AMR.h"
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
    // Length to dimensionalize nondim data
    amrex::Real dimL{0.};
    // Timestep size of reference data
    amrex::Real dt_modes{0.0};

    // Time (in HOS sim) to begin reading from
    // (this corresponds to t = 0 in AMR-Wind sim, functions as offset)
    amrex::Real t_winit{-1.0};
    // Timestep (in HOS sim) of latest imported wave data
    int ntime{0};
    // Timestep (in HOS sim) of first imported wave data
    int n_winit{0};
    // Timestep (in HOS sim) of last available wave data
    int n_wstop{0};
    // Timestep offset (in HOS sim) for resetting to earlier data
    int n_offset{0};

    // Time (in AMR-Wind sim) when ow arrays were last changed
    // (this is from the last timestep, goes with time-interpolated data)
    amrex::Real t_last{-1.0};
    // Time (in AMR-Wind sim) of latest imported wave data
    amrex::Real t{0.0};

    // Vectors to store modes
    std::vector<std::complex<double>> mX, mY, mZ, mFS;

// Struct variables that have special types
#ifdef AMR_WIND_USE_W2A
    // FFTW plan
    fftw_plan plan;
    // FFTW pointers to store modes ready for ifft
    fftw_complex *eta_mptr = nullptr, *u_mptr = nullptr, *v_mptr = nullptr,
                 *w_mptr = nullptr;
    // ReadModes object
    ReadModes rmodes;
#endif

    // Height vector (where velocity is sampled) stuff
    amrex::Vector<amrex::Real> hvec;
    // Index vector (where hvec overlaps with local boxes)
    amrex::Vector<int> indvec;
    // Flag indicating interpolation should take place on this processor
    bool do_interp{true};
    // Flag indicating regrid has occurred since the last resizing of vecs
    bool resize_flag{false};
    // Vectors for spatial data from transformed modes
    amrex::Gpu::DeviceVector<amrex::Real> sp_eta_vec, sp_u_vec, sp_v_vec,
        sp_w_vec;
    // Spacing of spatial data
    amrex::Real sp_dx{0.}, sp_dy{0.};

    ~W2AWavesData()
    {
#ifdef AMR_WIND_USE_W2A
        if (eta_mptr) {
            delete[] eta_mptr;
        }
        if (u_mptr) {
            delete[] u_mptr;
        }
        if (v_mptr) {
            delete[] v_mptr;
        }
        if (w_mptr) {
            delete[] w_mptr;
        }
        fftw_destroy_plan(plan);
#endif
    }
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
