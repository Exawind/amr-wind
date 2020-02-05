#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create the experimental setup
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_tuscan()
{
  /****************************************************************************
   * Define some variables                                                    *
   ***************************************************************************/

  Real xlo, xhi, ylo, yhi, zlo, zhi;

  Real xlen, ylen;
  Real xmid, ymid;

  Array<Real,3> point{0.0, 0.0, 0.0};
  Array<Real,3> normal{0.0, 0.0, 0.0};

  Real zlen = 0.2;

  xlo = 0.125;
  xhi = 0.875;

  ylo = xlo;
  yhi = xhi;

  zlo = zlen;
  zhi = 1.0 - zlen;

  xlen = xhi-xlo;
  ylen = yhi-ylo;

  xmid = xlo + xlen/2.0;
  ymid = ylo + ylen/2.0;

  /****************************************************************************
   * Build the lower and upper boxes                                          *
   ***************************************************************************/

  point  = {xmid,  ylo,  0.0};
  normal = { 0.0,  1.0,  0.0};
  EB2::PlaneIF b1(point, normal);

  point  = { xhi, ymid,  0.0};
  normal = {-1.0,  0.0,  0.0};
  EB2::PlaneIF b2(point, normal);

  point  = {xmid,  yhi,  0.0};
  normal = { 0.0, -1.0,  0.0};
  EB2::PlaneIF b3(point, normal);

  point  = { xlo, ymid,  0.0};
  normal = { 1.0,  0.0,  0.0};
  EB2::PlaneIF b4(point, normal);

  point  = {xmid, ymid,  zlo};
  normal = { 0.0,  0.0, -1.0};
  EB2::PlaneIF a1(point, normal);

  point  = {xmid, ymid,  zhi};
  normal = { 0.0,  0.0,  1.0};
  EB2::PlaneIF a2(point, normal);

  auto box1 = EB2::makeIntersection(b1, b2, b3, b4, a1);
  auto box2 = EB2::makeIntersection(b1, b2, b3, b4, a2);

  /****************************************************************************
   * Build the center box connector                                           *
   ***************************************************************************/

  Real middle_frac = 0.25;

  xlo = xlo + middle_frac*xlen;
  xhi = xhi - middle_frac*xlen;

  ylo = ylo + middle_frac*ylen;
  yhi = yhi - middle_frac*ylen;

  point  = {xmid,  ylo,  0.0};
  normal = { 0.0,  1.0,  0.0};
  EB2::PlaneIF c1(point, normal);

  point  = { xhi, ymid,  0.0};
  normal = {-1.0,  0.0,  0.0};
  EB2::PlaneIF c2(point, normal);

  point  = {xmid,  yhi,  0.0};
  normal = { 0.0, -1.0,  0.0};
  EB2::PlaneIF c3(point, normal);

  point  = { xlo, ymid,  0.0};
  normal = { 1.0,  0.0,  0.0};
  EB2::PlaneIF c4(point, normal);

  auto box3 = EB2::makeIntersection(c1, c2, c3, c4);

  // Merge all the boxes together.

  auto boxes = EB2::makeUnion(box1, box2, box3);

  //   /****************************************************************************
  //    *                                                                          *
  //    * Build standard EB Factories                                              *
  //    *                                                                          *
  //    ***************************************************************************/

  auto boxes_inside = EB2::makeComplement(boxes);
  auto gshop        = EB2::makeShop(boxes_inside);

  // Build index space
  int max_level_here = 0;
  int max_coarsening_level = 100;
  EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
