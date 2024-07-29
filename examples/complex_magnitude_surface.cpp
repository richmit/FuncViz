// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      complex_magnitude_surface.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Complex magnitude surface with coloring.@EOL
 @std       C++23
 @copyright 
  @parblock
  Copyright (c) 2024, Mitchell Jay Richling <http://www.mitchr.me/> All rights reserved.

  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software
     without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.
  @endparblock
 @filedetails   

  One popular way to plot complex functions is to use a surface plot of @f$\vert f(z)\vert@f$ and color the surface with @f$\arg(f(z))@f$.  This way we can
  simultaneously represent the magnitude and phase over the complex plain.  There are several ways to color the plots, and we will be following the method
  described by Richardson in 1991.  In this example, we demonstrates several techniques:

   - Alternate ways to do an initial sample (grid vs recursive)
   - Sample near a point in the domain
   - Sample below a level in the range
   - Sample near level curves
   - Sample based on a data value that is *not* part of the geometry
   - Sample near domain axis
   - Directly attach colors to geometric points
   - Do rough clipping with high sampling and cell filtering.  

  Eventually we will also demonstrate:

   - Clip with a clipping plain (TBD) -- requires new functionality in MR_rt_to_cc.
   - Extract level curves (TBD) -- requires new functionality in MR_rt_to_cc.

  References:
    Richardson (1991); Visualizing quantum scattering on the CM-2 supercomputer; Computer Physics Communications 63; pp 84-94"
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d9rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t cpf(tt_t::drpt_t xvec) {
  std::complex<double> z(xvec[0], xvec[1]);
  double z_abs, z_arg, f_re, f_im, f_abs, f_arg, red, green, blue;

  z_abs = std::abs(z);
  z_arg = std::arg(z);

  if ( (std::abs(z-1.0) > 1.0e-5) && (std::abs(z+1.0) > 1.0e-5) ) {
    std::complex<double> f;
    double f_abs2, f_re_scl, f_im_scl, f_abs2p1, ofs;
    f        = 1.0/(z+1.0) + 1.0/(z-1.0);
    f_re     = std::real(f);
    f_im     = std::imag(f);
    f_abs    = std::abs(f);
    f_arg    = std::arg(f);
    f_abs2   = f_abs * f_abs;
    f_re_scl = f_re / std::sqrt(30.0/5.0);
    f_im_scl = f_im / std::sqrt(2.0);
    f_abs2p1 = 1 + f_abs2;
    ofs      = (f_abs<1 ? -1.0 : 1.0) * (0.5 - f_abs/f_abs2p1);
    red      = ofs + (0.5 + (std::sqrt(2.0/3.0) * f_re) / f_abs2p1);
    green    = ofs + (0.5 - (f_re_scl - f_im_scl)       / f_abs2p1);
    blue     = ofs + (0.5 - (f_re_scl + f_im_scl)       / f_abs2p1);
  } else {
    f_re = f_im = f_abs = f_arg = red = green = blue = std::numeric_limits<double>::quiet_NaN();
  }

  return {z_abs, z_arg, f_re, f_im, f_abs, f_arg, red, green, blue};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::src_t  cpfd(tt_t::drpt_t xvec) {
  int    idx_for_z = 4;
  double cut_for_z = 3.5;
  auto   fv        = cpf(xvec);

  if(std::isnan(fv[idx_for_z]))
    return 100000.0;
  else
    return fv[idx_for_z]-cut_for_z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.2, -1.2}, 
            { 2.2,  1.2});
  cc_t ccplx;
  tc_t treeConverter;

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Initial sample

  // On a uniform grid
  tree.refine_grid(3, cpf);

  // Alternately we can use refine_recursive() instead (refine_grid() is faster)
  // tree.refine_recursive(4, cpf);

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Sample near 0+0i because we have a minimum at that piont

  // The most direct method 
  // tree.refine_leaves_recursive_cell_pred(6, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_close_to_domain_point({0, 0}, 1.0e-2, i)); });

  // This function is positive with a universal minimum at 0+0i, so we could just sample where  |f| is below 1/4
  tree.refine_leaves_recursive_cell_pred(6, cpf, [&tree](tt_t::diti_t i) { return !(tree.cell_above_range_level(i, 4, 0.25, 1.0e-5)); });

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Sample around the poles where we will clip the graph

  // With nice ranges the singularities will be precicely located on cell vertexes.  So we can just refine NaNs.
  // tree.refine_recursive_if_cell_vertex_is_nan(6, cpf);

  // Or we can directly sample on the clip level at |f|=3.5.  
  tree.refine_leaves_recursive_cell_pred(7, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_range_level(i, 4, 3.5)); });

  // We can do the above with a constructed SDF instead.
  // tree.refine_leaves_recursive_cell_pred(6, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, cpfd)); });

  // Just like the previous, but with atomic refinement.
  // tree.refine_leaves_atomically_if_cell_pred(6, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, cpfd)); });

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Refine where we plan to draw level curves

  // The easiest thing is to use cell_cross_range_level() for this.
  for(auto lev: {0.4, 0.7, 1.1, 1.4, 1.8, 2.6, 3.5}) 
    tree.refine_leaves_recursive_cell_pred(7, cpf, [&tree, lev](tt_t::diti_t i) { return (tree.cell_cross_range_level(i, 4, lev)); });

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // We will be coloring based on arg(f), and so want to sample near the abrubpt change near arg(f)=0.

  // We can do this just like the level curves with |f|, but use arg(f) instead -- i.e. index 5 instead of 4.
  tree.refine_leaves_recursive_cell_pred(7, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_range_level(i, 5, 0.0)); });

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // We can sample near the real & imagaxes axes.

  // Sample near the real axis
  tree.refine_leaves_recursive_cell_pred(5, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_domain_level(i, 0, 0.0, 1.0e-6)); });

  // Sample near the imaginary axis
  tree.refine_leaves_recursive_cell_pred(5, cpf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_domain_level(i, 1, 0.0, 1.0e-6)); });

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  // We don't need to balance the three, but it makes things look nice.

  // Balance the three to the traditional level of 1 (no  cell borders a cell more than half it's size)
  tree.balance_tree(1, cpf);

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  tree.dump_tree(5);

  auto tcret = treeConverter.construct_geometry_fans(ccplx,
                                                     tree,
                                                     tree.get_leaf_cells_pred(tree.ccc_get_top_cell(), 
                                                                              [&tree](tt_t::diti_t i) { return !(tree.cell_above_range_level(i, 4, 3.5, 1.0e-6)); }),
                                                     2,
                                                     {{tc_t::tree_val_src_t::DOMAIN, 0}, 
                                                      {tc_t::tree_val_src_t::DOMAIN, 1},
                                                      {tc_t::tree_val_src_t::RANGE,  4}});

  // Note the first argument need not name *every* data element, just the first ones.
  ccplx.create_named_datasets({"Re(z)", "Im(z)", "abs(z)", "arg(z)", "Re(f(z))", "Im(f(z))", "abs(f(z))", "arg(f(z))"}, {{"COLORS", {8, 9, 10}}});

  std::cout << "TC Return: " << tcret << std::endl;

  ccplx.dump_cplx(5);

  ccplx.write_legacy_vtk("complex_magnitude_surface.vtk", "complex_magnitude_surface");
  ccplx.write_xml_vtk(   "complex_magnitude_surface.vtu", "complex_magnitude_surface");
  ccplx.write_ply(       "complex_magnitude_surface.ply", "complex_magnitude_surface");
}
/** @endcond */
