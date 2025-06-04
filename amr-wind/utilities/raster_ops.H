#ifndef RASTER_OPS_H
#define RASTER_OPS_H

#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>
#include <string>
#include <fstream>
#include <sstream>

/**
 * Reference: ArcGIS ASCII raster format specification
 * https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm
 */

namespace amr_wind::utils {

/** \brief Minimal raster reader for ArcGIS ASCII raster files (.asc) on a
 * regular grid.
 *
 *  Stores only the lower-left corner, cell size, and value array.
 *  Provides a static, GPU-friendly, optimized bilinear interpolation.
 *
 *  Usage:
 *    RasterASC raster;
 *    raster.read("file.asc");
 *    // Host or device interpolation:
 *    amrex::Real val = RasterASC::interp_value(
 *        raster.x0(), raster.y0(), raster.dx(),
 *        raster.nx(), raster.ny(), raster.value_ptr(), x, y,
 * raster.nodata_value());
 *
 *  See:
 * https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm
 */
class RasterASC
{
public:
    RasterASC() = default;

    /**
     * \brief Optimized bilinear interpolation for regular grids (host/device).
     *
     * \param xll Lower-left x coordinate
     * \param yll Lower-left y coordinate
     * \param dx  Cell size (square, as per ArcGIS ASCII spec)
     * \param nx  Number of x grid points
     * \param ny  Number of y grid points
     * \param vals Pointer to flattened value array (row-major: [j * nx + i])
     * \param x   Query x coordinate
     * \param y   Query y coordinate
     * \param nodata The no data value
     * \return    Interpolated value at (x, y)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real interp_value(
        amrex::Real xll,
        amrex::Real yll,
        amrex::Real dx,
        int nx,
        int ny,
        const amrex::Real* vals,
        amrex::Real x,
        amrex::Real y,
        amrex::Real nodata)
    {
        amrex::Real fx = (x - xll) / dx;
        amrex::Real fy = (y - yll) / dx;
        int i = static_cast<int>(fx);
        int j = static_cast<int>(fy);

        // Clamp indices to valid range
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= nx - 1) i = nx - 2;
        if (j >= ny - 1) j = ny - 2;

        amrex::Real ddx = fx - i;
        amrex::Real ddy = fy - j;

        int idx00 = j * nx + i;
        int idx10 = j * nx + (i + 1);
        int idx01 = (j + 1) * nx + i;
        int idx11 = (j + 1) * nx + (i + 1);

        amrex::Real v00 = vals[idx00];
        amrex::Real v10 = vals[idx10];
        amrex::Real v01 = vals[idx01];
        amrex::Real v11 = vals[idx11];

        // Check for nodata at any corner
        if (v00 == nodata || v10 == nodata || v01 == nodata || v11 == nodata) {
            return nodata;
        }

        amrex::Real v0 = v00 * (1.0 - ddx) + v10 * ddx;
        amrex::Real v1 = v01 * (1.0 - ddx) + v11 * ddx;
        return v0 * (1.0 - ddy) + v1 * ddy;
    }

    /**
     * \brief Interpolate the raster value at (x, y) using the member data.
     *
     * \param x Query x coordinate
     * \param y Query y coordinate
     * \return  Interpolated value at (x, y)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    interp(amrex::Real x, amrex::Real y) const
    {
        return interp_value(
            m_xllcorner, m_yllcorner, m_dx, m_ncols, m_nrows, m_values.data(),
            x, y, m_nodata);
    }

    /// \brief Pointer to value array
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const amrex::Real*
    value_ptr() const
    {
        return m_values.data();
    }
    /// \brief Number of x grid points
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int nx() const { return m_ncols; }

    /// \brief Number of y grid points
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int ny() const { return m_nrows; }

    /// \brief Lower-left x coordinate
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real x0() const
    {
        return m_xllcorner;
    }

    /// \brief Lower-left y coordinate
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real y0() const
    {
        return m_yllcorner;
    }

    /// \brief Cell size (always square for ArcGIS ASCII)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real dx() const
    {
        return m_dx;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nodata_value() const
    {
        return m_nodata;
    };

    /**
     * \brief Read ArcGIS ASCII raster file (.asc) and populate value array.
     *
     * The file is expected to have a 6-line header followed by ncols*nrows
     * values. The values are stored in row-major order, starting from the top
     * row.
     *
     * Throws amrex::Abort on failure.
     *
     * \param filename Path to the .asc file
     */
    void read(const std::string& filename)
    {
        std::ifstream file(filename, std::ios::in);
        if (!file.good()) {
            amrex::Abort("Cannot open raster file: " + filename);
        }

        std::string key;
        if (!(file >> key >> m_ncols)) {
            amrex::Abort("Failed to read ncols");
        }
        if (!(file >> key >> m_nrows)) {
            amrex::Abort("Failed to read nrows");
        }
        if (!(file >> key >> m_xllcorner)) {
            amrex::Abort("Failed to read xllcorner");
        }
        if (!(file >> key >> m_yllcorner)) {
            amrex::Abort("Failed to read yllcorner");
        }
        if (!(file >> key >> m_dx)) {
            amrex::Abort("Failed to read cellsize");
        }
        if (!(file >> key >> m_nodata)) {
            amrex::Abort("Failed to read nodata");
        }

        m_values.resize(m_ncols * m_nrows);

        // ArcGIS ASCII starts from top row, so we read from last y to first
        int count = 0;
        for (int j = m_nrows - 1; j >= 0; --j) {
            for (int i = 0; i < m_ncols; ++i) {
                if (!(file >> m_values[j * m_ncols + i])) {
                    amrex::Abort(
                        "Failed to read raster value at (" + std::to_string(i) +
                        "," + std::to_string(j) + ")");
                }
                ++count;
            }
        }
        if (count != m_ncols * m_nrows) {
            amrex::Abort(
                "Raster file does not contain the expected number of values");
        }
        if (!file.good() && !file.eof()) {
            amrex::Abort("Error occurred while reading raster file");
        }
    }

    /// Test-only: set all member data directly (for unit tests)
    void from_data(
        int ncols,
        int nrows,
        amrex::Real xll,
        amrex::Real yll,
        amrex::Real dx,
        amrex::Real nodata,
        const amrex::Vector<amrex::Real>& values)
    {
        m_ncols = ncols;
        m_nrows = nrows;
        m_xllcorner = xll;
        m_yllcorner = yll;
        m_dx = dx;
        m_nodata = nodata;
        m_values = values;
    }

private:
    amrex::Vector<amrex::Real> m_values;
    int m_ncols{0}, m_nrows{0};
    amrex::Real m_xllcorner{0.0}, m_yllcorner{0.0};
    amrex::Real m_dx{1.0};
    amrex::Real m_nodata{-9999.0};
};

} // namespace amr_wind::utils

#endif /* RASTER_OPS_H */