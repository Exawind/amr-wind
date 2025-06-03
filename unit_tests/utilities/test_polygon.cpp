#include "amr-wind/utilities/tagging/Polygon.H"
#include "gtest/gtest.h"

namespace amr_wind_tests {

using amr_wind::polygon_utils::Polygon;

TEST(Polygon, BasicInclusion)
{
    Polygon poly;
    // Square: (0,0) -> (0,10) -> (10,10) -> (10,0) -> (0,0)
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 0.0});
    poly.compute_bounding_box();

    EXPECT_TRUE(poly.contains({5.0, 5.0}));
    EXPECT_FALSE(poly.contains({15.0, 5.0}));
    EXPECT_FALSE(poly.contains({0.0, 5.0})); // on edge
}

TEST(Polygon, WithHole)
{
    Polygon poly;
    // Outer square
    poly.add_outer_vertex({0.0, 0.0});
    poly.add_outer_vertex({0.0, 10.0});
    poly.add_outer_vertex({10.0, 10.0});
    poly.add_outer_vertex({10.0, 0.0});
    poly.add_outer_vertex({0.0, 0.0});
    // Hole: (2,2) -> (2,8) -> (8,8) -> (8,2) -> (2,2)
    poly.start_hole();
    poly.add_vertex({2.0, 2.0});
    poly.add_vertex({2.0, 8.0});
    poly.add_vertex({8.0, 8.0});
    poly.add_vertex({8.0, 2.0});
    poly.add_vertex({2.0, 2.0});
    poly.compute_bounding_box();

    EXPECT_FALSE(poly.contains({5.0, 5.0})); // in hole
    EXPECT_TRUE(poly.contains({1.0, 1.0}));  // in shell
}

TEST(Polygon, Empty)
{
    Polygon poly;
    EXPECT_TRUE(poly.is_empty());
    EXPECT_FALSE(poly.contains({0.0, 0.0}));
}

} // namespace amr_wind_tests