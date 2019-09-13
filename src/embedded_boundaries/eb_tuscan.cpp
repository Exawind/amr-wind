#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

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

  xlo = 0.125;
  xhi = 0.875;

  ylo = xlo;
  yhi = xhi;

  zlo = 0.3125;
  zhi = 0.6875;

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

  xlo = xlo + 0.25*xlen;
  xhi = xhi - 0.25*xlen;

  ylo = ylo + 0.25*ylen;
  yhi = yhi - 0.25*ylen;

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

    auto gshop = EB2::makeShop(boxes);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EBSupport m_eb_support_level = EBSupport::full;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();

    // Make the EBFabFactory
    for(int lev = 0; lev <= max_level; lev++)
    {
        const EB2::Level& eb_is_lev = eb_is.getLevel(geom[lev]);
        eb_level = &eb_is_lev;
        ebfactory[lev].reset(new EBFArrayBoxFactory(*eb_level,
                                                    geom[lev],
                                                    grids[lev],
                                                    dmap[lev],
                                                    {m_eb_basic_grow_cells,
                                                    m_eb_volume_grow_cells,
                                                    m_eb_full_grow_cells},
                                                    m_eb_support_level));
    }
}
