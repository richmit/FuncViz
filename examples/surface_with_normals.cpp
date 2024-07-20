// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_with_normals.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Simple surface function graph.@EOL
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

  This example demonstrates several things
   - Surface normals
     - How to compute a surface gradient for a function plot
     - How to unitize the gradient into a surface normal
     - How to add the normal to the sample data stored by a MRPTree
   - Signed Distance Functions
     - How to create a simple SDF function
     - How to use an SDF to increase sampling on the SDF boundary
   - How to balance a tree
   - Tree conversion into geometry
     - How to convert a tree of sample data into a cell complex
     - How to check the result of treeConverter
   - How to dump a cell complex into various file types
*/
/*******************************************************************************************************************************************************.H.E.**/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::array<double, 4> dampCosWave2(std::array<double, 2> xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double d = x*x+y*y;
  double m = std::exp(-d/4);
  double s = std::sqrt(d);
  double z = m*cos(4*s);
  double dx = -(cos((4 * s)) * s + 4 * sin( (4 * s))) * x * exp(-x * x / 2 - y * y / 2);
  double dy = -(cos((4 * s)) * s + 4 * sin( (4 * s))) * y * exp(-x * x / 2 - y * y / 2);
  if (s>1.0e-5) {
    dx = dx / s;
    dy = dy / s;
  } else {
    dx = 0;
    dy = 0;
  }
  double nm = std::sqrt(1+dx*dx+dy*dy);
  return {z, -dx/nm, -dy/nm, 1/nm};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double circleSDF2(double r, std::array<double, 2> xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double m = x*x+y*y;
  return (r*r-m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d4rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.1, -2.1}, 
            { 2.1,  2.1});
  cc_t ccplx;
  tc_t treeConverter;

  // Make a few samples on a uniform grid
  tree.refine_grid(2, dampCosWave2);

  // Sample more on the humps (we could have done this with a derivative test, but this is a good example of how to use an SDF)
  for(double i: {0, 1, 2, 3}) {
    double r = i*std::numbers::pi/4;
    tree.refine_leaves_recursive_cell_pred(7, dampCosWave2, [&tree, r](int i) { return (tree.cell_cross_sdf(i, std::bind_front(circleSDF2, r))); });
  }

  // Balance the three to the traditional level of 1 (no  cell borders a cell more than half it's size)
  // tree.balance_tree(1, dampCosWave2);

  tree.dump_tree(5);

  treeConverter.construct_geometry(ccplx,
                                   tree,
                                   tc_t::cell_structure_t::FANS, 
                                   2,
                                   { "points", 
                                     tc_t::tree_val_src_t::DOMAIN, 0, 
                                     tc_t::tree_val_src_t::DOMAIN, 1,
                                     tc_t::tree_val_src_t::RANGE,  0},
                                   {{ "x",      tc_t::tree_val_src_t::DOMAIN, 0},
                                    { "y",      tc_t::tree_val_src_t::DOMAIN, 1},
                                    { "f(x,y)", tc_t::tree_val_src_t::RANGE,  0},
                                    { "-df/dx", tc_t::tree_val_src_t::RANGE,  1},
                                    { "-df/dy", tc_t::tree_val_src_t::RANGE,  2}},
                                   {{"NORMALS",
                                     tc_t::tree_val_src_t::RANGE, 1,
                                     tc_t::tree_val_src_t::RANGE, 2,
                                     tc_t::tree_val_src_t::RANGE, 3}});

  ccplx.dump_cplx(5);
  ccplx.write_legacy_vtk("surface_with_normals.vtk", "surface_with_normals");
  ccplx.write_xml_vtk(   "surface_with_normals.vtu", "surface_with_normals");
  ccplx.write_ply(       "surface_with_normals.ply", "surface_with_normals");
}
