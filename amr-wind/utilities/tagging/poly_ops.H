#ifndef POLY_OPS_H
#define POLY_OPS_H

#include "AMReX.H"
#include <cmath>
#include "AMReX_REAL.H"
#include "AMReX_Gpu.H"

// Taken from https://www.geeksforgeeks.org/point-in-polygon-in-cpp/
namespace amr_wind::polygon_utils {

struct Point
{
    amrex::Real x, y;
};

// Function to compute the cross product of vectors (p1p2)
// and (p1p3)

AMREX_FORCE_INLINE amrex::Real
cross_product(const Point& p1, const Point& p2, const Point& p3)
{
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

// Function to check if point p lies on segment p1p2
AMREX_FORCE_INLINE bool
is_point_on_segment(const Point& p, const Point& p1, const Point& p2)
{
    // Check if point p lies on the line segment p1p2 and
    // within the bounding box of p1p2
    return cross_product(p1, p2, p) == 0 && p.x >= std::min(p1.x, p2.x) &&
           p.x <= std::max(p1.x, p2.x) && p.y >= std::min(p1.y, p2.y) &&
           p.y <= std::max(p1.y, p2.y);
}

// Function to compute the winding number of a point with
// respect to a polygon
int winding_number(const amrex::Vector<Point>& polygon, const Point& point)
{
    int n = polygon.size();
    int windingNumber = 0;

    // Iterate through each edge of the polygon
    for (int i = 0; i < n; i++) {
        Point p1 = polygon[i];
        // Next vertex in the polygon
        Point p2 = polygon[(i + 1) % n];

        // Check if the point lies on the current edge
        if (is_point_on_segment(point, p1, p2)) {
            // Point is on the polygon boundary
            return 0;
        }

        // Calculate the cross product to determine winding
        // direction
        if (p1.y <= point.y) {
            if (p2.y > point.y && cross_product(p1, p2, point) > 0) {
                windingNumber++;
            }
        } else {
            if (p2.y <= point.y && cross_product(p1, p2, point) < 0) {
                windingNumber--;
            }
        }
    }
    // Return the winding number
    return windingNumber;
}

// Function to check if a point is inside a polygon using
// the winding number algorithm
bool is_point_in_polygon(
    const amrex::Vector<Point>& polygon, const Point& point)
{
    // Compute the winding number for the point with respect
    // to the polygon
    return winding_number(polygon, point) != 0;
}

} // namespace amr_wind::polygon_utils

#endif /* POLY_OPS_H */
