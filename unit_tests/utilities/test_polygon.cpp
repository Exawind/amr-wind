#include "amr-wind/utilities/tagging/Polygon.H"
#include <gtest/gtest.h>
#include <iostream>

using namespace amr_wind::polygon_utils;

TEST(Polygon, SimpleSquare)
{
    Polygon poly;
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 0.0});
    poly.compute_bounding_box();

    const auto* pts = poly.ring_ptr(0);
    int n_outer = poly.ring_size(0);
    const auto* all_pts = poly.points().data();
    const auto* ring_offsets = poly.ring_offsets().data();
    int n_rings = poly.num_rings();
    int n_points = poly.num_points();

    std::cout << "SimpleSquare: Testing point (5,5) (inside)\n";
    EXPECT_TRUE(poly.contains({5.0, 5.0}));
    EXPECT_TRUE(Polygon::is_point_in_ring(pts, n_outer, {5.0, 5.0}));
    EXPECT_TRUE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {5.0, 5.0}));

    std::cout << "SimpleSquare: Testing point (0,5) (on edge)\n";
    EXPECT_FALSE(poly.contains({0.0, 5.0}));
    EXPECT_FALSE(Polygon::is_point_in_ring(pts, n_outer, {0.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {0.0, 5.0}));

    std::cout << "SimpleSquare: Testing point (-1,5) (outside)\n";
    EXPECT_FALSE(poly.contains({-1.0, 5.0}));
    EXPECT_FALSE(Polygon::is_point_in_ring(pts, n_outer, {-1.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {-1.0, 5.0}));
}

TEST(Polygon, PolygonWithHole)
{
    Polygon poly;
    // Outer square
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 0.0});
    poly.start_hole();
    // Hole: smaller square
    poly.add_vertex({3.0, 3.0});
    poly.add_vertex({3.0, 7.0});
    poly.add_vertex({7.0, 7.0});
    poly.add_vertex({7.0, 3.0});
    poly.add_vertex({3.0, 3.0});
    poly.compute_bounding_box();

    const auto* all_pts = poly.points().data();
    const auto* ring_offsets = poly.ring_offsets().data();
    int n_rings = poly.num_rings();
    int n_points = poly.num_points();

    std::cout << "PolygonWithHole: Testing point (2,2) (inside outer, outside "
                 "hole)\n";
    EXPECT_TRUE(poly.contains({2.0, 2.0}));
    EXPECT_TRUE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {2.0, 2.0}));

    std::cout << "PolygonWithHole: Testing point (5,5) (inside hole)\n";
    EXPECT_FALSE(poly.contains({5.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {5.0, 5.0}));

    std::cout << "PolygonWithHole: Testing point (10,5) (on exterior edge)\n";
    EXPECT_FALSE(poly.contains({10.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {10.0, 5.0}));

    std::cout << "PolygonWithHole: Testing point (3,5) (on hole edge)\n";
    EXPECT_FALSE(poly.contains({3.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {3.0, 5.0}));

    std::cout << "PolygonWithHole: Testing point (-1,5) (outside)\n";
    EXPECT_FALSE(poly.contains({-1.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {-1.0, 5.0}));
}

TEST(Polygon, ComplexPolygonSelfIntersecting)
{
    Polygon poly;
    // Bowtie (self-intersecting)
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({5.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({0.0, 0.0});
    poly.compute_bounding_box();

    const auto* all_pts = poly.points().data();
    const auto* ring_offsets = poly.ring_offsets().data();
    int n_rings = poly.num_rings();
    int n_points = poly.num_points();

    std::cout << "ComplexPolygonSelfIntersecting: Testing point (5,5) (center, "
                 "ambiguous)\n";
    EXPECT_FALSE(poly.contains({5.0, 5.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {5.0, 5.0}));
}

TEST(Polygon, DegenerateCases)
{
    Polygon poly;
    // Single point
    poly.add_outer_vertex({1.0, 1.0});
    poly.compute_bounding_box();
    const auto* all_pts = poly.points().data();
    const auto* ring_offsets = poly.ring_offsets().data();
    int n_rings = poly.num_rings();
    int n_points = poly.num_points();

    std::cout << "DegenerateCases: Testing single point polygon at (1,1)\n";
    EXPECT_FALSE(poly.contains({1.0, 1.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {1.0, 1.0}));

    // Line
    Polygon poly2;
    poly2.add_outer_vertex({0.0, 0.0});
    poly2.add_outer_vertex({1.0, 1.0});
    poly2.compute_bounding_box();

    const auto* all_pts2 = poly2.points().data();
    const auto* ring_offsets2 = poly2.ring_offsets().data();
    int n_rings2 = poly2.num_rings();
    int n_points2 = poly2.num_points();

    std::cout << "DegenerateCases: Testing line polygon from (0,0) to (1,1)\n";
    EXPECT_FALSE(poly2.contains({0.5, 0.5}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts2, ring_offsets2, n_rings2, n_points2, {0.5, 0.5}));
}

TEST(Polygon, BoundingBox)
{
    Polygon poly;
    poly.add_outer_vertex({-2.0, 3.0});
    poly.add_outer_vertex({4.0, 5.0});
    poly.add_outer_vertex({1.0, -1.0});
    poly.add_outer_vertex({-2.0, 3.0});
    poly.compute_bounding_box();

    Polygon::Point bbox_lo = {0.0, 0.0};
    Polygon::Point bbox_hi = {0.0, 0.0};
    poly.get_bounding_box(bbox_lo, bbox_hi);

    std::cout << "BoundingBox: Testing point (0,0) (inside bbox)\n";
    EXPECT_TRUE(poly.bounding_box_contains({0.0, 0.0}));
    EXPECT_TRUE(
        Polygon::poly_bounding_box_contains({0.0, 0.0}, bbox_lo, bbox_hi));
    std::cout << "BoundingBox: Testing point (10,10) (outside bbox)\n";
    EXPECT_FALSE(poly.bounding_box_contains({10.0, 10.0}));
    EXPECT_FALSE(
        Polygon::poly_bounding_box_contains({10.0, 10.0}, bbox_lo, bbox_hi));
}

TEST(Polygon, MultipleHoles)
{
    Polygon poly;
    // Outer square
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 0.0});
    poly.start_hole();
    // Hole 1
    poly.add_vertex({2.0, 2.0});
    poly.add_vertex({2.0, 4.0});
    poly.add_vertex({4.0, 4.0});
    poly.add_vertex({4.0, 2.0});
    poly.add_vertex({2.0, 2.0});
    poly.start_hole();
    // Hole 2
    poly.add_vertex({6.0, 6.0});
    poly.add_vertex({6.0, 8.0});
    poly.add_vertex({8.0, 8.0});
    poly.add_vertex({8.0, 6.0});
    poly.add_vertex({6.0, 6.0});
    poly.compute_bounding_box();

    const auto* all_pts = poly.points().data();
    const auto* ring_offsets = poly.ring_offsets().data();
    int n_rings = poly.num_rings();
    int n_points = poly.num_points();

    std::cout
        << "MultipleHoles: Testing point (5,5) (inside outer, outside holes)\n";
    EXPECT_TRUE(poly.contains({5.0, 5.0}));
    EXPECT_TRUE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {5.0, 5.0}));

    std::cout << "MultipleHoles: Testing point (3,3) (inside hole 1)\n";
    EXPECT_FALSE(poly.contains({3.0, 3.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {3.0, 3.0}));

    std::cout << "MultipleHoles: Testing point (7,7) (inside hole 2)\n";
    EXPECT_FALSE(poly.contains({7.0, 7.0}));
    EXPECT_FALSE(
        Polygon::is_point_in_polygon(
            all_pts, ring_offsets, n_rings, n_points, {7.0, 7.0}));
}