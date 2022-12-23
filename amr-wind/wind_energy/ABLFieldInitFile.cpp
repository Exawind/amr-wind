#include <cmath>

#include "amr-wind/wind_energy/ABLFieldInitFile.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

namespace amr_wind {

ABLFieldInitFile::ABLFieldInitFile()
{
    amrex::ParmParse pp_abl("ABL");
    // Get netcdf input file name
    pp_abl.get("initial_condition_input_file", m_ic_input);
}

bool ABLFieldInitFile::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const int lev) const
{
    // Initialized using lev to prevent compiler warnings when netcdf is not
    // used. Will always initialize to false. cppcheck-suppress constVariable
    bool interp_fine_levels = (lev < 0);

#ifdef AMR_WIND_USE_NETCDF
    // Skip fine levels and interpolate data from already loaded coarse levels
    if (!m_ic_input.empty() && lev > 0) {
        interp_fine_levels = true;
    }
    // Load the netcdf file with data if specified in the inputs
    if (!m_ic_input.empty() && lev == 0) {

        // Open the netcdf input file
        // This file should have the same dimensions as the simulation
        auto ncf = ncutils::NCFile::open(m_ic_input, NC_NOWRITE);

        // Ensure that the input dimensions match the coarsest grid size
        const auto& domain = geom.Domain();

        // The indices that determine the start and end points of the i, j, k
        // arrays. The max and min are there to ensure that the points form
        // ghost cells are not used
        auto i0 = amrex::max(vbx.smallEnd(0), domain.smallEnd(0));
        auto i1 = amrex::min(vbx.bigEnd(0), domain.bigEnd(0));

        auto j0 = amrex::max(vbx.smallEnd(1), domain.smallEnd(1));
        auto j1 = amrex::min(vbx.bigEnd(1), domain.bigEnd(1));

        auto k0 = amrex::max(vbx.smallEnd(2), domain.smallEnd(2));
        auto k1 = amrex::min(vbx.bigEnd(2), domain.bigEnd(2));

        // The x, y and z velocity components (u, v, w)
        auto uvel = ncf.var("uvel");
        auto vvel = ncf.var("vvel");
        auto wvel = ncf.var("wvel");

        // Loop through all points in the domain and set velocities to values
        // from the input file
        // start is the first index from where to read data
        amrex::Vector<size_t> start{
            {static_cast<size_t>(i0), static_cast<size_t>(j0),
             static_cast<size_t>(k0)}};
        // count is the total number of elements to read in each direction
        amrex::Vector<size_t> count{
            {static_cast<size_t>(i1 - i0 + 1), static_cast<size_t>(j1 - j0 + 1),
             static_cast<size_t>(k1 - k0 + 1)}};

        // Working vector to read data onto host
        amrex::Vector<double> tmp;
        int dlen = count[0] * count[1] * count[2];
        tmp.resize(dlen);
        // Vector to store the 3d data into a single array and set size
        amrex::Gpu::DeviceVector<amrex::Real> uvel_d(dlen, 0.0);
        amrex::Gpu::DeviceVector<amrex::Real> vvel_d(dlen, 0.0);
        amrex::Gpu::DeviceVector<amrex::Real> wvel_d(dlen, 0.0);

        // Read the velocity components u, v, w and copy to device
        uvel.get(tmp.data(), start, count);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, &tmp[0], &tmp[dlen - 1] + 1,
            uvel_d.begin());
        vvel.get(tmp.data(), start, count);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, &tmp[0], &tmp[dlen - 1] + 1,
            vvel_d.begin());
        wvel.get(tmp.data(), start, count);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, &tmp[0], &tmp[dlen - 1] + 1,
            wvel_d.begin());

        // Pointers to velocity objects
        const auto* uvel_dptr = uvel_d.data();
        const auto* vvel_dptr = vvel_d.data();
        const auto* wvel_dptr = wvel_d.data();

        // Get count components for device
        int ct1 = count[1];
        int ct2 = count[2];

        // Amrex parallel for to assign the velocity at each point
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // The counter to go from 3d to 1d vector
                auto idx = (i - i0) * ct2 * ct1 + (j - j0) * ct1 + (k - k0);
                // Pass values from temporary array to the velocity field
                velocity(i, j, k, 0) = uvel_dptr[idx];
                velocity(i, j, k, 1) = vvel_dptr[idx];
                velocity(i, j, k, 2) = wvel_dptr[idx];
            });
        // Close the netcdf file
        ncf.close();
    }
#endif
    return interp_fine_levels;
}

} // namespace amr_wind
