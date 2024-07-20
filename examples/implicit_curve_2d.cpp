// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      implicit_curve_2d.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Sampleing on a 2D grid to extract an implicit curve.@EOL
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

  Many visualization tools can extract a "level set" from a triangulation.  In this example we extract a curve from a function of two variables -- usually
  called an "implicit curve".  The trick to obtaining high quality results is to make sure the triangulation has a high enough resolution.  The naive way, and
  the way it's normally done, is to simply sample uniformly across a grid the curve is expected to be inside.  A better way is to detect where the curve is,
  and to sample at higher resolution near the curve.
*/
/*******************************************************************************************************************************************************.H.E.**/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is a classic "difficult case" for implicit curve algorithms.
double f(std::array<double, 2> xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double z = ((2*x*x*y - 2*x*x - 3*x + y*y*y - 33*y + 32) * ((x-2)*(x-2) + y*y + 3))/3000;
  if (z>1.0)
    z=1.0;
  if (z<-1.0)
    z=-1.0;
  return z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d1rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-10.0, -6.5},
            { 10.0,  6.5});
  cc_t ccplx;
  tc_t treeConverter;

  // First we sample the top cell.  Just one cell!
  tree.sample_cell(f);

  // Now we recursively refine cells that seem to cross over the curve
  tree.refine_leaves_recursive_cell_pred(7, f, [&tree](tt_t::diti_t i) { return (tree.cell_cross_range_level(i, 0, 0.0)); });

  // We could have used the function f as an SDF, and achieved the same result with the following:
  // tree.refine_leaves_recursive_cell_pred(7, f, [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, f)); });

  tree.dump_tree(20);

  // Convert the geometry into a 3D dataset so we can see the contour on the surface
  treeConverter.construct_geometry(ccplx,
                                   tree,
                                   tc_t::cell_structure_t::FANS, 
                                   2,
                                   { "points", 
                                     tc_t::tree_val_src_t::DOMAIN, 0,
                                     tc_t::tree_val_src_t::DOMAIN, 1,
                                     tc_t::tree_val_src_t::RANGE,  0},
                                   {{"x",       tc_t::tree_val_src_t::DOMAIN, 0},
                                    {"y",       tc_t::tree_val_src_t::DOMAIN, 1},
                                    {"f(x,y)",  tc_t::tree_val_src_t::RANGE,  0}},
                                   {});

  ccplx.write_xml_vtk("implicit_curve_2d.vtu", "implicit_curve_2d");
}

