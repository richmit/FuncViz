// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      geomi_seg_isect_type.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Unit tests for MR_cell_cplx.@EOL
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
*/
/*******************************************************************************************************************************************************.H.E.**/

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>

#include "MR_cell_cplx.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(geomi_seg_isect_type) {

  mjr::MRccT5 aPoly;

  aPoly.add_point({0.0,  0.0, 0.0}); // 0: 
  aPoly.add_point({5.0,  0.0, 0.0}); // 1: 0-1=[0,5]
  aPoly.add_point({6.0,  0.0, 0.0}); // 2: 1-2=[5,6]  0-2=[0,6]
  aPoly.add_point({8.0,  0.0, 0.0}); // 3: 0-3=[0,8]  2-3=[6,8]

  aPoly.add_point({0.0, -1.0, 0.0}); // 4: /-
  aPoly.add_point({5.0,  1.0, 0.0}); // 5: /+

  aPoly.add_point({6.0, -1.0, 0.0}); // 6: /-
  aPoly.add_point({7.0,  1.0, 0.0}); // 7: /+
  aPoly.add_point({8.0,  2.0, 0.0}); // 8: /++

  // ---------------------------------------------------------------------------
  //                                      0  1  2  3  4  5  6  7  8
  // ---------------------------------------------------------------------------
  // 0,0,0,1 => [0,0] [0,5] BAD_SEGMENT   *
  //                                      [              ]
  //                                      *
  // ---------------------------------------------------------------------------
  // 0,1,2,3 => [0,5] [6,8] C0_EMPTY      [              ]
  //                                                        [     ] 
  // ---------------------------------------------------------------------------
  // 0,1,1,2 => [0,5] [5,6] C1_VERTEX1    [              ]
  //                                                     [  ] 
  //                                                     *
  // ---------------------------------------------------------------------------
  // 0,3,4,5 => [0,8] [-,+] C1_INTERIOR   [                       ] 
  //                                      /              /
  //                                             *
  // ---------------------------------------------------------------------------
  // 0,1,0,1 => [0,5] [0,5] CI_VERTEX2    [              ]
  //                                      [              ] 
  //                                      ****************
  // ---------------------------------------------------------------------------
  // 0,1,0,2 => [0,5] [0,6] CI_VERTEX1    [              ]
  //                                      ****************
  //                                      [                 ] 
  // ---------------------------------------------------------------------------
  // 0,3,1,2 => [0,8] [5,6] CI_VERTEX0    [                       ]
  //                                                     [  ] 
  //                                                     ****
  // ---------------------------------------------------------------------------
  // 0,2,1,3 => [0,6] [5,8] CI_VERTEX0    [                 ]
  //                                                     [        ] 
  //                                                     ****
  // ---------------------------------------------------------------------------
  // 0,3,0,5 => [0,8] [0,+] C1_VERTEX1    [                       ] 
  //                                      [              /
  //                                      *
  // ---------------------------------------------------------------------------
  // 0,1,2,5 => [0,5] [6,+] C0_EMPTY      [              ] 
  //                                                        [              /
  // ---------------------------------------------------------------------------
  // 0,1,6,7 => [0,5] [+,-] C0_EMPTY      [              ]  
  //                                                     /  /
  // ---------------------------------------------------------------------------
  // 0,1,6,7 => [0,5] [+,-] C0_EMPTY      [              ]  
  //                                                           /  /
  // ---------------------------------------------------------------------------
  // 0,1,5,7 => [0,5] [+,-] C0_EMPTY      [              ]  
  //                                                     <     >
  // ---------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------

  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 0, 0, 1) == mjr::MRccT5::seg_isect_t::BAD_SEGMENT);
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 2, 3) == mjr::MRccT5::seg_isect_t::C0_EMPTY);        // colinear case
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 2, 5) == mjr::MRccT5::seg_isect_t::C0_EMPTY);        // non-colinear (1 pt in second segment colinear)
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 6, 7) == mjr::MRccT5::seg_isect_t::C0_EMPTY);        // non-colinear (no end points of second seg colinear with first, sec seg straddles line containing first segment)
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 7, 8) == mjr::MRccT5::seg_isect_t::C0_EMPTY);        // non-colinear (no end points of second seg colinear with first, sec seg dosen't straddle line containing first segment)
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 5, 7) == mjr::MRccT5::seg_isect_t::C0_EMPTY);        // non-colinear (segments parallel)
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 1, 2) == mjr::MRccT5::seg_isect_t::C1_VERTEX1);      // colinear case
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 3, 0, 5) == mjr::MRccT5::seg_isect_t::C1_VERTEX1);      // non-colinear
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 3, 4, 5) == mjr::MRccT5::seg_isect_t::C1_INTERIOR);
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 0, 1) == mjr::MRccT5::seg_isect_t::CI_VERTEX2);
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 1, 0, 2) == mjr::MRccT5::seg_isect_t::CI_VERTEX1); 
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 3, 1, 2) == mjr::MRccT5::seg_isect_t::CI_VERTEX0);
  BOOST_CHECK(aPoly.geomi_seg_isect_type(0, 2, 1, 3) == mjr::MRccT5::seg_isect_t::CI_VERTEX0);

}
