// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      vector_field_3d.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Vector field for Lorenz system.@EOL
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

  This example illustrates how to uniformly sample a vector field.  Just for fun we have also produced a solution to the Lorenz system, and directly
  stored it with a MR_cell_cplx.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::array<double, 3> vf(std::array<double, 3> xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double z = xvec[2];
  double a = 10.0;
  double b = 28.0;
  double c = 8.0/3.0;
  return { a*y-a*z,
           x*b-x*z,
           x*y-c*z
         };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b3d3rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t vftree({-30.0, -30.0,  -0.0},
              { 30.0,  30.0,  60.0});
  cc_t vfccplx;
  tc_t vftreeConverter;

  /* Uniform sampling */
  vftree.refine_grid(5, vf);

  /* Dump the vector field */
  vftreeConverter.construct_geometry(vfccplx,
                                     vftree,
                                     tc_t::cell_structure_t::RECTANGLES, 
                                     0,
                                     { "points", 
                                       tc_t::tree_val_src_t::DOMAIN,  0,
                                       tc_t::tree_val_src_t::DOMAIN,  1,
                                       tc_t::tree_val_src_t::DOMAIN,  2},
                                     {{"x",   tc_t::tree_val_src_t::DOMAIN, 0},
                                      {"y",   tc_t::tree_val_src_t::DOMAIN, 1},
                                      {"z",   tc_t::tree_val_src_t::DOMAIN, 2}},
                                     {{"d",
                                       tc_t::tree_val_src_t::RANGE,  0,
                                       tc_t::tree_val_src_t::RANGE,  1,
                                       tc_t::tree_val_src_t::RANGE,  2}});
  vfccplx.dump_cplx(5);
  vfccplx.write_xml_vtk("vector_field_3d-f.vtu", "vector_field_3d-f");


  cc_t cvccplx;

  int max_steps = 100000;
  double delta  = 0.001;
  double t      = 0;
  double x_old  = 0.1;
  double y_old  = 0.0;
  double z_old  = 0.0;
  double a      = 10.0;
  double b      = 28.0;
  double c      = 8.0 / 3.0;

  auto p_old = cvccplx.add_point({x_old, y_old, z_old});
  cvccplx.add_data_if_new("t", t);
  for(int num_steps=0;num_steps<max_steps;num_steps++) {    
    double x_new = x_old + a*(y_old-x_old)*delta;
    double y_new = y_old + (x_old*(b-z_old)-y_old)*delta;
    double z_new = z_old + (x_old*y_old-c*z_old)*delta;
    t += delta;
    auto p_new = cvccplx.add_point({x_new, y_new, z_new});
    cvccplx.add_data_if_new("t", t);
    cvccplx.add_cell(cc_t::cell_type_t::SEGMENT, {p_old, p_new});
    x_old=x_new;
    y_old=y_new;
    z_old=z_new;
    p_old=p_new;

  }

  cvccplx.dump_cplx(5);
  cvccplx.write_xml_vtk("vector_field_3d-c.vtu", "vector_field_3d-c");


}
/** @endcond */
