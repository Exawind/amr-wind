/**
 * \file raster_ops.H
 * \brief Minimal raster reader and interpolator for ArcGIS ASCII raster files
 * (.asc) on a regular grid.
 *
 * Reference: ArcGIS ASCII raster format specification
 * https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm
 *
 * This class supports both cell-centered ("area") and node-centered ("point")
 * pixel interpretations. By default, the raster is interpreted as cell-centered
 * ("area").
 *
 * Features:
 *  - Read ArcGIS ASCII raster files (.asc)
 *  - Set raster data directly for testing
 *  - Bilinear interpolation and nearest-neighbor lookup
 *  - Supports both cell-centered and node-centered grids
 *  - Handles nodata values
 *
 * Usage Example:
 * \code
 * #include "amr-wind/utilities/raster_ops.H"
 * using amr_wind::utils::RasterASC;
 * using amr_wind::utils::PixelInterpretation;
 * using amr_wind::utils::pixel_interpretation_from_string;
 *
 * // --- Host usage ---
 * RasterASC raster;
 * raster.read("myraster.asc");
 * amrex::Real val = raster.interp(10.5, 20.5);
 * amrex::Real nearest = raster.find_nearest(10.5, 20.5);
 *
 * // --- Static methods (host or device/GPU) ---
 * // Suppose you have a device pointer to values and want to interpolate on the
 * GPU: const amrex::Real* vals = raster.value_ptr(); int nx = raster.nx(); int
 * ny = raster.ny(); amrex::Real xll = raster.x0(); amrex::Real yll =
 * raster.y0(); amrex::Real dx = raster.dx(); amrex::Real nodata =
 * raster.nodata_value(); PixelInterpretation interp =
 * raster.pixel_interpretation();
 *
 * // Example: Interpolate at many points on the GPU
 * amrex::Vector<amrex::Real> xs = ...; // coordinates on device
 * amrex::Vector<amrex::Real> ys = ...;
 * amrex::Vector<amrex::Real> out_vals(xs.size());
 *
 * amrex::ParallelFor(xs.size(), [=] AMREX_GPU_DEVICE (int i) noexcept {
 *     out_vals[i] = RasterASC::interp_value(
 *         xll, yll, dx, nx, ny, vals, xs[i], ys[i], nodata, interp);
 * });
 *
 * // Example: Find nearest value at many points on the GPU
 * amrex::ParallelFor(xs.size(), [=] AMREX_GPU_DEVICE (int i) noexcept {
 *     out_vals[i] = RasterASC::find_nearest_value(
 *         xll, yll, dx, nx, ny, vals, xs[i], ys[i], nodata, interp);
 * });
 *
 * // --- Setting pixel interpretation ---
 * raster.set_pixel_interpretation(PixelInterpretation::Point);
 * // or using a string:
 * raster.set_pixel_interpretation(pixel_interpretation_from_string("area"));
 *
 * // --- Set data directly (for tests or synthetic rasters) ---
 * amrex::Vector<amrex::Real> vals_vec = {1.0, 2.0, 3.0, 4.0};
 * raster.from_data(2, 2, 0.0, 0.0, 1.0, -9999.0, vals_vec,
 * PixelInterpretation::Area);
 * \endcode
 *
 * PixelInterpretation options:
 *   - Area: cell-centered (default, Arc/ESRI ASCII "area")
 *   - Point: node-centered
 *
 * To convert from string: use pixel_interpretation_from_string("area") or
 * ("point").
 */

#ifndef RASTER_OPS_H
#define RASTER_OPS_H

#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>
#include <string>
#include <fstream>
#include <sstream>

namespace amr_wind::utils {

enum class PixelInterpretation {
    Area, // cell-centered (default, Arc/ESRI ASCII "area")
    Point // node-centered
};

inline PixelInterpretation
pixel_interpretation_from_string(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "area" || lower == "cell" || lower == "cellcentered" ||
        lower == "cell-centered") {
        return PixelInterpretation::Area;
    }
    if (lower == "point" || lower == "node" || lower == "nodecentered" ||
        lower == "node-centered") {
        return PixelInterpretation::Point;
    }
    amrex::Abort("Unknown PixelInterpretation string: " + s);
    return PixelInterpretation::Area; // Silence warning
}

inline std::string to_string(PixelInterpretation interp)
{
    switch (interp) {
    case PixelInterpretation::Area:
        return "area";
    case PixelInterpretation::Point:
        return "point";
    default:
        return "unknown";
    }
}

class RasterASC
{
public:
    RasterASC() = default;

    void set_pixel_interpretation(PixelInterpretation interp)
    {
        m_pixel_interp = interp;
    }

    /**
     * \brief Find the value at the nearest grid point to (x, y) (host/device).
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
     * \param pixel_interp Pixel interpretation (Area = cell-centered, Point =
     * node-centered)
     * \return    Value at nearest grid point to (x, y)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    find_nearest_value(
        amrex::Real xll,
        amrex::Real yll,
        amrex::Real dx,
        int nx,
        int ny,
        const amrex::Real* vals,
        amrex::Real x,
        amrex::Real y,
        amrex::Real nodata,
        PixelInterpretation pixel_interp = PixelInterpretation::Area)
    {
        // Compute indices (assuming (xll, yll) is the lower-left corner, and
        // (0,0) is the first cell)
        amrex::Real fx = (x - xll) / dx;
        amrex::Real fy = (y - yll) / dx;
        if (pixel_interp == PixelInterpretation::Area) {
            fx -= 0.5;
            fy -= 0.5;
        }
        int i = static_cast<int>(amrex::Math::round(fx));
        int j = static_cast<int>(amrex::Math::round(fy));

        // Clamp indices to valid range
        if (i < 0) {
            i = 0;
        }
        if (j < 0) {
            j = 0;
        }
        if (i >= nx) {
            i = nx - 1;
        }
        if (j >= ny) {
            j = ny - 1;
        }

        amrex::Real val = vals[j * nx + i];
        if (val == nodata) {
            return nodata;
        }
        return val;
    }

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
     * \param pixel_interp Pixel interpretation (Area = cell-centered, Point =
     * node-centered)
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
        amrex::Real nodata,
        PixelInterpretation pixel_interp = PixelInterpretation::Area)
    {
        amrex::Real fx = (x - xll) / dx;
        amrex::Real fy = (y - yll) / dx;
        if (pixel_interp == PixelInterpretation::Area) {
            fx -= 0.5;
            fy -= 0.5;
        }
        int i = static_cast<int>(fx);
        int j = static_cast<int>(fy);

        // Clamp indices to valid range
        if (i < 0) {
            i = 0;
        }
        if (j < 0) {
            j = 0;
        }
        if (i >= nx - 1) {
            i = nx - 2;
        }
        if (j >= ny - 1) {
            j = ny - 2;
        }

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
            x, y, m_nodata, m_pixel_interp);
    }

    /**
     * \brief Find the value at the nearest grid point to (x, y)
     *
     * \param x Query x coordinate
     * \param y Query y coordinate
     * \return  Value at nearest grid point to (x, y)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    find_nearest(amrex::Real x, amrex::Real y) const
    {
        return find_nearest_value(
            m_xllcorner, m_yllcorner, m_dx, m_ncols, m_nrows, m_values.data(),
            x, y, m_nodata, m_pixel_interp);
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

    /// \brief Get the nodata value
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nodata_value() const
    {
        return m_nodata;
    };

    /// \brief Get the current pixel interpretation
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PixelInterpretation
    pixel_interpretation() const
    {
        return m_pixel_interp;
    }

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

        m_values.resize(
            static_cast<std::size_t>(m_ncols) *
            static_cast<std::size_t>(m_nrows));

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

    /**
     * \brief Set all raster metadata and values directly (for tests or
     * synthetic rasters).
     *
     * \param ncols         Number of columns (x grid points)
     * \param nrows         Number of rows (y grid points)
     * \param xll           Lower-left x coordinate
     * \param yll           Lower-left y coordinate
     * \param dx            Cell size (square, as per ArcGIS ASCII spec)
     * \param nodata        The no data value
     * \param values        Flattened value array (row-major: [j * ncols + i])
     * \param interp        Pixel interpretation (Area = cell-centered, Point =
     * node-centered)
     */
    void from_data(
        int ncols,
        int nrows,
        amrex::Real xll,
        amrex::Real yll,
        amrex::Real dx,
        amrex::Real nodata,
        const amrex::Vector<amrex::Real>& values,
        PixelInterpretation interp = PixelInterpretation::Area)
    {
        m_ncols = ncols;
        m_nrows = nrows;
        m_xllcorner = xll;
        m_yllcorner = yll;
        m_dx = dx;
        m_nodata = nodata;
        m_values = values;
        m_pixel_interp = interp;
    }

private:
    PixelInterpretation m_pixel_interp{PixelInterpretation::Area};
    amrex::Vector<amrex::Real> m_values;
    int m_ncols{0}, m_nrows{0};
    amrex::Real m_xllcorner{0.0}, m_yllcorner{0.0};
    amrex::Real m_dx{1.0};
    amrex::Real m_nodata{-9999.0};
};

} // namespace amr_wind::utils

#endif /* RASTER_OPS_H */