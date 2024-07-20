// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      performance_with_large_surface.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Stress test with a large surface object.@EOL
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

  Just a nice parametric surface without any weirdness.  Because we use a large mesh for this one, I have including timing code.  Note c(u,v) can be used to
  render stripes on the surface.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <chrono>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::array<double, 4> shellStripes2(std::array<double, 2> xvec) {
  double u = std::numbers::pi   * xvec[0] + std::numbers::pi + 0.1;
  double v = std::numbers::pi/2 * xvec[1] + std::numbers::pi/2;
  return {u*sin(u)*cos(v), 
          u*std::cos(u)*std::cos(v), 
          u*sin(v), 
          std::fmod(u*sin(v), 2)
         };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d4rT            tt_t;
typedef mjr::MRccT5                  cc_t;   // Replace with mjr::MRccF5, and compare treeConverter performance.
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::system_clock::now();
  tt_t tree;
  cc_t ccplx;
  tc_t treeConverter;
  std::chrono::time_point<std::chrono::system_clock> constructTime = std::chrono::system_clock::now();

  tree.refine_grid(9, shellStripes2);
  std::chrono::time_point<std::chrono::system_clock> sampleTime = std::chrono::system_clock::now();

  tree.dump_tree(20);
  std::chrono::time_point<std::chrono::system_clock> dumpTime = std::chrono::system_clock::now();

  treeConverter.construct_geometry(ccplx,
                                   tree,
                                   tc_t::cell_structure_t::FANS, 
                                   2,
                                   { "points", 
                                     tc_t::tree_val_src_t::RANGE,  0,
                                     tc_t::tree_val_src_t::RANGE,  1,
                                     tc_t::tree_val_src_t::RANGE,   2},
                                   {{"u",       tc_t::tree_val_src_t::DOMAIN, 0},
                                    {"v",       tc_t::tree_val_src_t::DOMAIN, 1},
                                    {"x(u,v)",  tc_t::tree_val_src_t::RANGE,  0},
                                    {"y(u,v)",  tc_t::tree_val_src_t::RANGE,  1},
                                    {"z(u,v)",  tc_t::tree_val_src_t::RANGE,  2},
                                    {"c(u,v)",  tc_t::tree_val_src_t::RANGE,  3}},
                                   {});
  std::chrono::time_point<std::chrono::system_clock> vtkFanTime = std::chrono::system_clock::now();

  ccplx.write_xml_vtk("performance_with_large_surface.vtu", "performance_with_large_surface");
  std::chrono::time_point<std::chrono::system_clock> vtkWriteTime = std::chrono::system_clock::now();

  std::cout << "constructTime time ... " << static_cast<std::chrono::duration<double>>(constructTime-startTime)     << " sec" << std::endl;
  std::cout << "sampleTime time ...... " << static_cast<std::chrono::duration<double>>(sampleTime-constructTime)    << " sec" << std::endl;
  std::cout << "dumpTime time ........ " << static_cast<std::chrono::duration<double>>(dumpTime-sampleTime)         << " sec" << std::endl;
  std::cout << "treeConverter time ... " << static_cast<std::chrono::duration<double>>(vtkFanTime-dumpTime)         << " sec" << std::endl;
  std::cout << "write_vtk time ....... " << static_cast<std::chrono::duration<double>>(vtkWriteTime-vtkFanTime)     << " sec" << std::endl;
  std::cout << "Total Run Time ....... " << static_cast<std::chrono::duration<double>>(vtkWriteTime-startTime)      << " sec" << std::endl;
}
/** @endcond */
