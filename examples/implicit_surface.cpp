// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      implicit_surface.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Sampling for an implicit surface.@EOL
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

  This example is very similar to implicit_curve_2d.cpp; however, instead of extracting a curve from a triangulation of a surface, this time we extract a
  surface from a quad tessellation of a hexahedron.
*/
/*******************************************************************************************************************************************************.H.E.**/

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
