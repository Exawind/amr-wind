#ifndef AMR_WIND_POLYGON_H_
#define AMR_WIND_POLYGON_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <limits>
#include <ostream>

namespace amr_wind::polygon_utils {

/**
 * \class Polygon
 * \ingroup amr_utils
 * \brief 2D polygon with support for holes and point-in-polygon tests.
 *
 * Stores all points in a single array, with offsets for each ring (outer and
 * holes). This layout is GPU-friendly and matches usage in TerrainRefinement.
 *
 * Example input for a polygon with one hole:
 * \code
 * tagging.myrefine.type = TerrainRefinement
 * tagging.myrefine.vertical_distance = 100
 * tagging.myrefine.level = 1
 * tagging.myrefine.poly_exterior = 0 0  0 10  10 10  10 0  0 0
 * tagging.myrefine.poly_num_holes = 1
 * tagging.myrefine.poly_hole_0 = 2 2  2 8  8 8  8 2  2 2
 * \endcode
 *
 * - poly_exterior: List of x y pairs for the outer ring (must be closed).
 * - poly_num_holes: Number of holes.
 * - poly_hole_N: List of x y pairs for each hole (must be closed).
 */
class Polygon
{
public:
    using Point = amrex::Array<amrex::Real, 2>;

    Polygon() = default;

    /** \brief Returns true if the polygon has no points */
    bool is_empty() const { return m_points.empty(); }

    /**
     * \brief Get the bounding box of the polygon
     *
     * \param[out] bbox_min Lower left corner of bounding box (x, y)
     * \param[out] bbox_max Upper right corner of bounding box (x, y)
     */
    void get_bounding_box(Point& bbox_min, Point& bbox_max) const
    {
        bbox_min = m_bbox_min;
        bbox_max = m_bbox_max;
    }

    /**
     * \brief Returns true if the bounding box has been computed and is valid
     *
     * The bounding box is considered valid if min < max in both x and y.
     * This is typically set by calling compute_bounding_box().
     *
     * \return True if bounding box is valid, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool bbox_valid() const
    {
        return (m_bbox_min[0] < m_bbox_max[0]) &&
               (m_bbox_min[1] < m_bbox_max[1]);
    }

    /** \brief Add a vertex to the outer ring (must be called before any holes)
     */
    void add_outer_vertex(const Point& pt)
    {
        if (m_ring_offsets.empty()) {
            m_ring_offsets.push_back(0);
        }
        m_points.push_back(pt);
    }

    /** \brief Start a new hole (inner ring) */
    void start_hole()
    {
        m_ring_offsets.push_back(static_cast<int>(m_points.size()));
    }

    /** \brief Add a vertex to the current ring (outer or last hole) */
    void add_vertex(const Point& pt) { m_points.push_back(pt); }

    /** \brief Return the number of rings (outer + holes) */
    int num_rings() const { return static_cast<int>(m_ring_offsets.size()); }

    /** \brief Return the number of points in ring i */
    int ring_size(int i) const
    {
        if (i < 0 || i >= static_cast<int>(m_ring_offsets.size())) {
            return 0;
        }
        int start = m_ring_offsets[i];
        int end = (i + 1 < static_cast<int>(m_ring_offsets.size()))
                      ? m_ring_offsets[i + 1]
                      : static_cast<int>(m_points.size());
        return end - start;
    }

    /** \brief Return a pointer to the start of ring i */
    const Point* ring_ptr(int i) const
    {
        if (i < 0 || i >= static_cast<int>(m_ring_offsets.size())) {
            return nullptr;
        }
        return m_points.data() + m_ring_offsets[i];
    }

    /**
     * \brief Check if a point is inside the polygon (excluding holes)
     * \param pt The point to test
     * \return True if inside, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    contains(const Point& pt) const
    {
        if (is_empty()) {
            return false;
        }
        if (!bounding_box_contains(pt)) {
            return false;
        }
        return is_point_in_polygon(
            m_points.data(), m_ring_offsets.data(), num_rings(), num_points(),
            pt);
    }

    /**
     * \brief Standalone check if a point is inside a bounding box
     * \param pt The point to test
     * \param bbox_min Lower left corner of bounding box
     * \param bbox_max Upper right corner of bounding box
     * \return True if inside bounding box, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
    poly_bounding_box_contains(
        const Point& pt, const Point& bbox_min, const Point& bbox_max)
    {
        return (
            pt[0] >= bbox_min[0] && pt[0] <= bbox_max[0] &&
            pt[1] >= bbox_min[1] && pt[1] <= bbox_max[1]);
    }

    /**
     * \brief Check if a point is inside the bounding box
     * \param pt The point to test
     * \return True if inside bounding box, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    bounding_box_contains(const Point& pt) const
    {
        return poly_bounding_box_contains(pt, m_bbox_min, m_bbox_max);
    }

    /** \brief Compute the bounding box of the polygon */
    void compute_bounding_box()
    {
        m_bbox_min = {
            std::numeric_limits<amrex::Real>::max(),
            std::numeric_limits<amrex::Real>::max()};
        m_bbox_max = {
            std::numeric_limits<amrex::Real>::lowest(),
            std::numeric_limits<amrex::Real>::lowest()};
        for (const auto& pt : m_points) {
            m_bbox_min[0] = amrex::min(m_bbox_min[0], pt[0]);
            m_bbox_min[1] = amrex::min(m_bbox_min[1], pt[1]);
            m_bbox_max[0] = amrex::max(m_bbox_max[0], pt[0]);
            m_bbox_max[1] = amrex::max(m_bbox_max[1], pt[1]);
        }
    }

    /**
     * \brief Read polygon data from ParmParse
     * \param prefix ParmParse prefix for polygon input
     */
    void read_from_parmparse(const std::string& prefix)
    {
        amrex::ParmParse pp(prefix);

        // Outer ring
        amrex::Vector<amrex::Real> outer_coords;
        pp.queryarr("poly_exterior", outer_coords);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            outer_coords.size() % 2 == 0,
            "Polygon exterior ring must have an even number of coordinates.");

        m_points.clear();
        m_ring_offsets.clear();
        m_ring_offsets.push_back(0);
        for (amrex::Long i = 0;
             i < static_cast<amrex::Long>(outer_coords.size()); i += 2) {
            m_points.push_back({outer_coords[i], outer_coords[i + 1]});
        }

        // Holes
        int num_holes = 0;
        pp.query("poly_num_holes", num_holes);
        for (int h = 0; h < num_holes; ++h) {
            std::string key = "poly_hole_" + std::to_string(h);
            amrex::Vector<amrex::Real> hole_coords;
            pp.queryarr(key.c_str(), hole_coords);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                hole_coords.size() % 2 == 0,
                "Polygon hole must have an even number of coordinates.");

            m_ring_offsets.push_back(static_cast<int>(m_points.size()));
            for (amrex::Long i = 0;
                 i < static_cast<amrex::Long>(hole_coords.size()); i += 2) {
                m_points.push_back({hole_coords[i], hole_coords[i + 1]});
            }
        }

        compute_bounding_box();
    }

    /**
     * \brief Print the polygon to the given stream
     * \param os Output stream (defaults to amrex::OutStream())
     */
    void print(std::ostream& os = amrex::OutStream()) const
    {
        if (is_empty()) {
            os << "Polygon: (empty)\n";
            os << "Bounding box: not computed\n";
            return;
        }
        os << "Polygon (outer ring):\n";
        for (int i = 0; i < ring_size(0); ++i) {
            const auto& p = ring_ptr(0)[i];
            os << "  (" << p[0] << ", " << p[1] << ")\n";
        }
        for (int h = 1; h < num_rings(); ++h) {
            os << "  Hole " << h << ":\n";
            for (int i = 0; i < ring_size(h); ++i) {
                const auto& p = ring_ptr(h)[i];
                os << "    (" << p[0] << ", " << p[1] << ")\n";
            }
        }

        bool bbox_valid =
            (m_bbox_min[0] < m_bbox_max[0]) && (m_bbox_min[1] < m_bbox_max[1]);
        os << "Bounding box: ";
        if (bbox_valid) {
            os << "min=(" << m_bbox_min[0] << ", " << m_bbox_min[1] << "), "
               << "max=(" << m_bbox_max[0] << ", " << m_bbox_max[1] << ")\n";
        } else {
            os << "not computed\n";
        }
    }

    /** \brief Access to raw points for GPU */
    const amrex::Vector<Point>& points() const { return m_points; }

    /** \brief Get total number of points */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int num_points() const
    {
        return static_cast<int>(m_points.size());
    }

    /** \brief Access to ring offsets for GPU */
    const amrex::Vector<int>& ring_offsets() const { return m_ring_offsets; }
    /**
     * \brief Standalone ring point-in-polygon test for GPU use
     * \param ring Pointer to ring points
     * \param n Number of points in ring
     * \param pt Point to test
     * \return True if inside, false otherwise
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
    is_point_in_ring(const Point* ring, int n, const Point& pt)
    {
        int wn = 0;
        for (int j = 0; j < n; ++j) {
            const auto& p1 = ring[j];
            const auto& p2 = ring[(j + 1) % n];
            amrex::Real cross = (p2[0] - p1[0]) * (pt[1] - p1[1]) -
                                (pt[0] - p1[0]) * (p2[1] - p1[1]);
            if (amrex::Math::abs(cross) < 1e-12 &&
                pt[0] >= amrex::min(p1[0], p2[0]) &&
                pt[0] <= amrex::max(p1[0], p2[0]) &&
                pt[1] >= amrex::min(p1[1], p2[1]) &&
                pt[1] <= amrex::max(p1[1], p2[1])) {
                return false; // on edge = not inside
            }
            if (p1[1] <= pt[1]) {
                if (p2[1] > pt[1] && cross > 0) {
                    ++wn;
                }
            } else {
                if (p2[1] <= pt[1] && cross < 0) {
                    --wn;
                }
            }
        }
        return wn != 0;
    }

    /**
     * \brief Standalone check if a point is inside a polygon with holes, given
     * raw data pointers
     * \param p_poly_points Pointer to all polygon points
     * \param p_ring_offsets Pointer to ring offsets
     * \param n_rings Number of rings (outer + holes)
     * \param n_points Total number of points
     * \param pt The point to test
     * \return True if strictly inside, false otherwise (false on edge or in
     * hole)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool is_point_in_polygon(
        const Point* p_poly_points,
        const int* p_ring_offsets,
        int n_rings,
        int n_points,
        const Point& pt)
    {
        // Outer ring
        int start = static_cast<int>(p_ring_offsets[0]);
        int end = (n_rings > 1) ? static_cast<int>(p_ring_offsets[1])
                                : static_cast<int>(n_points);
        int n = end - start;
        // Check if pt is on edge of outer ring
        for (int j = 0; j < n; ++j) {
            const auto& p1 = p_poly_points[start + j];
            const auto& p2 = p_poly_points[start + ((j + 1) % n)];
            if (is_point_on_segment(pt, p1, p2)) {
                return false;
            }
        }
        if (!is_point_in_ring(p_poly_points + start, n, pt)) {
            return false;
        }
        // Check holes
        for (int ring_i = 1; ring_i < n_rings; ++ring_i) {
            int h_start = static_cast<int>(p_ring_offsets[ring_i]);
            int h_end = (ring_i + 1 < n_rings)
                            ? static_cast<int>(p_ring_offsets[ring_i + 1])
                            : static_cast<int>(n_points);
            int h_n = h_end - h_start;
            // Check if pt is on edge of hole
            for (int j = 0; j < h_n; ++j) {
                const auto& p1 = p_poly_points[h_start + j];
                const auto& p2 = p_poly_points[h_start + ((j + 1) % h_n)];
                if (is_point_on_segment(pt, p1, p2)) {
                    return false;
                }
            }
            if (is_point_in_ring(p_poly_points + h_start, h_n, pt)) {
                return false;
            }
        }
        return true;
    }

private:
    amrex::Vector<Point> m_points;
    amrex::Vector<int> m_ring_offsets;
    Point m_bbox_min;
    Point m_bbox_max;

    /** \brief Helper: check if point pt is on segment p1-p2 */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static bool
    is_point_on_segment(const Point& pt, const Point& p1, const Point& p2)
    {
        amrex::Real cross = is_left(p1, p2, pt);
        if (std::abs(cross) > 1e-12) {
            return false;
        }
        return (
            pt[0] >= amrex::min(p1[0], p2[0]) &&
            pt[0] <= amrex::max(p1[0], p2[0]) &&
            pt[1] >= amrex::min(p1[1], p2[1]) &&
            pt[1] <= amrex::max(p1[1], p2[1]));
    }

    /** \brief Compute winding number for a point with respect to ring i */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
    winding_number(int ring_idx, const Point& pt) const
    {
        if (ring_idx < 0 ||
            ring_idx >= static_cast<int>(m_ring_offsets.size())) {
            return 0;
        }
        int wn = 0;
        int start = m_ring_offsets[ring_idx];
        int end = (ring_idx + 1 < static_cast<int>(m_ring_offsets.size()))
                      ? m_ring_offsets[ring_idx + 1]
                      : static_cast<int>(m_points.size());
        int n = end - start;
        if (n == 0) {
            return 0;
        }
        for (int j = 0; j < n; ++j) {
            const auto& p1 = m_points[start + j];
            const auto& p2 = m_points[start + ((j + 1) % n)];
            if (is_point_on_segment(pt, p1, p2)) {
                return 0;
            }
            if (p1[1] <= pt[1]) {
                if (p2[1] > pt[1] && is_left(p1, p2, pt) > 0) {
                    ++wn;
                }
            } else {
                if (p2[1] <= pt[1] && is_left(p1, p2, pt) < 0) {
                    --wn;
                }
            }
        }
        return wn;
    }

    /** \brief Helper for winding number: is pt left of line p0->p1? */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    is_left(const Point& p0, const Point& p1, const Point& p2)
    {
        return (p1[0] - p0[0]) * (p2[1] - p0[1]) -
               (p2[0] - p0[0]) * (p1[1] - p0[1]);
    }
};

} // namespace amr_wind::polygon_utils

#endif // AMR_WIND_POLYGON_H_