// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      parametric_curve_3d.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Parametric curve as the intersection of two parametric surfaces.@EOL
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

  This example illustrates a simple parametric curve in 3D.  Stuff illistrated:
   - Produce a VTK file representing a parametric curve in 3D.
   - Produce a VTK file representing a parametric surface in 3D -- actually two VTK files each with a diffrent surface.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b1d3rT              tt1_t;
typedef mjr::MRccT5                    cc1_t;
typedef mjr::MR_rt_to_cc<tt1_t, cc1_t> tc1_t;

typedef mjr::tree15b2d3rT              tt2_t;
typedef mjr::MRccT5                    cc2_t;
typedef mjr::MR_rt_to_cc<tt2_t, cc2_t> tc2_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt1_t::rrpt_t twisted_cubic_crv(tt1_t::drpt_t t) {
  return { t, t*t, t*t*t };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::rrpt_t twisted_cubic_srf1(tt2_t::drpt_t uvvec) {
  tt2_t::src_t u = uvvec[0];
  tt2_t::src_t v = uvvec[1];
  return { u, u*u, v };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::rrpt_t twisted_cubic_srf2(tt2_t::drpt_t uvvec) {
  tt2_t::src_t u = uvvec[0];
  tt2_t::src_t v = uvvec[1];
  return { u, v, u*u*u };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt1_t crv_tree;
  cc1_t crv_ccplx;
  tc1_t crv_tree_conv;
  crv_tree.refine_grid(8, twisted_cubic_crv);
  crv_tree_conv.construct_geometry_fans(crv_ccplx,
                                        crv_tree,
                                        1,
                                        {{tc1_t::tree_val_src_t::RANGE, 0},
                                         {tc1_t::tree_val_src_t::RANGE, 1},
                                         {tc1_t::tree_val_src_t::RANGE, 2}});
  crv_ccplx.create_named_datasets({"t", "x(t)", "y(t)", "z(t)"});
  crv_ccplx.dump_cplx(5);
  crv_ccplx.write_xml_vtk("parametric_curve_3d-crv.vtu", "parametric_curve_3d-crv");

  tt2_t srf1_tree;
  cc2_t srf1_ccplx;
  tc2_t srf1_tree_conv;
  srf1_tree.refine_grid(5, twisted_cubic_srf1);
  srf1_tree_conv.construct_geometry_fans(srf1_ccplx,
                                         srf1_tree,
                                         2,
                                         {{tc2_t::tree_val_src_t::RANGE, 0},
                                          {tc2_t::tree_val_src_t::RANGE, 1},
                                          {tc2_t::tree_val_src_t::RANGE, 2}});
  srf1_ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  srf1_ccplx.dump_cplx(5);
  srf1_ccplx.write_xml_vtk("parametric_curve_3d-srf1.vtu", "parametric_curve_3d-srf1");

  tt2_t srf2_tree;
  cc2_t srf2_ccplx;
  tc2_t srf2_tree_conv;
  srf2_tree.refine_grid(5, twisted_cubic_srf2);
  srf2_tree_conv.construct_geometry_fans(srf2_ccplx,
                                         srf2_tree,
                                         2,
                                         {{tc2_t::tree_val_src_t::RANGE, 0},
                                          {tc2_t::tree_val_src_t::RANGE, 1},
                                          {tc2_t::tree_val_src_t::RANGE, 2}});
  srf2_ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  srf2_ccplx.dump_cplx(5);
  srf2_ccplx.write_xml_vtk("parametric_curve_3d-srf2.vtu", "parametric_curve_3d-srf2");
}
/** @endcond */
