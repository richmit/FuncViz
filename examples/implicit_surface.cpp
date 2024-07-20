////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double isf(std::array<double, 3> xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double z = xvec[2];
  return x*x*y+y*y*x-z*z*z-1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b3d1rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.3, -2.3, -2.3}, 
            { 2.3,  2.3,  2.3});
  cc_t ccplx;
  tc_t treeConverter;

  /* Initial uniform sample */
  tree.refine_grid(4, isf);

  /* Refine near surface */
  tree.refine_leaves_recursive_cell_pred(6, isf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, isf)); });

  /* Convert our tree to a cell complex */
  treeConverter.construct_geometry(ccplx,
                                   tree,
                                   tc_t::cell_structure_t::RECTANGLES, 
                                   3,
                                   { "points", 
                                     tc_t::tree_val_src_t::DOMAIN, 0, 
                                     tc_t::tree_val_src_t::DOMAIN, 1,
                                     tc_t::tree_val_src_t::DOMAIN, 2},
                                   {{ "x",        tc_t::tree_val_src_t::DOMAIN, 0},
                                    { "y",        tc_t::tree_val_src_t::DOMAIN, 1},
                                    { "z",        tc_t::tree_val_src_t::DOMAIN, 2},
                                    { "f(x,y,z)", tc_t::tree_val_src_t::RANGE,  0}},
                                   {});

  /* Display some data about the cell complex */
  ccplx.dump_cplx(5);

  /* Write out our cell complex */
  ccplx.write_xml_vtk("implicit_surface.vtu", "implicit_surface");
}
