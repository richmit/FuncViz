// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_cell_cplx.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Implementation for the MR_cell_cplx class.@EOL
 @keywords  VTK polydata PLY file MR_rect_tree polygon triangulation cell complex tessellation
 @std       C++23
 @see       MR_rect_tree.hpp, MR_rt_to_cc.hpp
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef MJR_INCLUDE_MR_cell_cplx

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>                                                     /* STL algorithm           C++11    */
#include <array>                                                         /* array template          C++11    */
#include <set>                                                           /* STL set                 C++98    */
#include <unordered_map>                                                 /* STL hash map            C++11    */
#include <map>                                                           /* STL map                 C++11    */
#include <vector>                                                        /* STL vector              C++11    */ 
#include <iostream>                                                      /* C++ iostream            C++11    */
#include <string>                                                        /* C++ strings             C++11    */
#include <fstream>                                                       /* C++ fstream             C++98    */
#include <iomanip>                                                       /* C++ stream formatting   C++11    */
#include <cmath>                                                         /* std:: C math.h          C++11    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put everything in the mjr namespace
namespace mjr {
  /** @brief Template Class used hold tessellations of MR_rect_tree geometric data as well as MR_rect_tree point data sets.

      This class's primary use case is to hold tessellations of MR_rect_tree geometric data as well as MR_rect_tree point data sets.  For many applications the
      final goal is to produce a file (usually an unstructured VTK file).  That said, this class is generally useful by itself.

      This code is a quick-n-dirty hack job.  I simply wanted to quickly get to a point where I could visually test MR_rect_tree objects.  Logically this
      class should be split into several classes (unique point list, unique cell list, 3D analytic geometry, etc...).  It also has some significant
      limitations:

        - VTK files generated from MR_rect_tree objects don't require all of the capability of VTK files
          - Only Unstructured_Grid files are supported
          - Only four types of cells are supported: points, segments, triangles, quads, pyramids, & hexahedrons.
            - No tetrahedrons!  This may seem bazaar; however, they are simply not a natural product of MR_rect_tree tessellation.
          - Only point data is supported (i.e. no cell data), and only scalar & vector data types are supported
          - Both legacy & XML files are supported, but they are ASCII only.
          - XML files are serial and self contained
          - NaN's in legacy files are not properly handled by many VTK applications
        - Performance isn't a priority
          - Especially true for 0-cells
          - Geometry checks are slow (may be avoided by using add_cell without cell_type_t argument)
          - Point de-duplication is slow (may be turned off via uniq_points template parameter)
          - Cell de-duplication is slow (may be turned off via uniq_cells template parameter)
          - In general expect a 5x speedup by turning off de-duplication & geometry checks.
        - Things we don't check
          - Memory allocation -- you run out, the thing crashes
          - Cells that are a part of other cells may be added (i.e. a segment that is part of an existing triangle)

      Cells supported:

       - Point
         - Zero-dimensional cell. 
         - Defined by a single point.
       
       - Segment
         - One-dimensional cell. 
         - Defined an ordered list of two points. 
         - The direction along the line is from the first point to the second point.
       
       - Triangle
         - Two-dimensional cell. 
         - Defined by an ordered list of three points. 
         - The points are ordered counterclockwise, defining a surface normal using the right-hand rule.
       
       - Quadrilateral
         - Two-dimensional cell. 
         - Defined by an ordered list of four points.
         - The points are ordered counterclockwise, defining a surface normal using the right-hand rule.
         - The points must be coplainer
         - The quadrilateral is convex and its edges must not intersect. 
       
       - Hexahedron
         - Three-dimensional cell consisting of six quadrilateral faces, twelve edges, and eight vertices. 
         - The hexahedron is defined by an ordered list of eight points. 
         - The faces and edges must not intersect any other faces and edges, and the hexahedron must be convex.
       
       - Pyramid
         - Three-dimensional cell consisting of one quadrilateral face, four triangular faces, eight edges, and five vertices. 
         - The pyramid is defined by an ordered list of five points. 
         - The four points defining the quadrilateral base plane must be convex
         - The fifth point, the apex, must not be co-planar with the base points.

         \verbatim                                                                                                              
                                                                                                                    4                 
                                                                                                                   .^.                  
                                                                              3---------2                       .  / \ .              
                                                                             /|        /|                     .   /   \  .                  
                                                                            / |       / |                  .     /     \   .          
                                           2            3---------2        7---------6  |                 4...../.......\.....0  Back 
                                          / \           |         |        |  |      |  |                 |    /         \    |       
          0        0-------1             /   \          |         |        |  0------|--1 Back            |   /           \   |       
                                        /     \         |         |        | /       | /                  |  /             \  |       
                                       0-------1        |         |        |/        |/                   | /               \ |       
                                                        0---------1        4---------5 Front              |/                 \|       
                                                                                                          2 ----------------- 3  Front
         \endverbatim

      A number of quality checks may be performed on points and cells before they are added to the object.  These checks can slow down execution of add_cell()
      by an order of magnitude, and almost double the RAM required for this class.  These checks are entirely optional.

        - uniq_points:
          - Many 3D file formats store the master list of potential points, and then use an integer index into this list when defining geometric objects.  
          - Most visualization software is pretty tolerant of having duplicate points on this list when just doing rendering.
          - Duplicate points can break many software packages pretty badly when it comes to computation.
          - Depending on how points are added to this object, this check can avoid quite a bit of wasted RAM and produce *much* smaller output files.
          - I normally leave this check turned on.
        - uniq_cells:
          - Most packages can render objects with duplicate cells, but you might see graphical artifacts or glitching.
          - When using a mesh with duplicate cells, computations can sometimes fail or take a *very* long time
          - I normally leave this check turned on.
        - chk_cell_vertexes:
          - This check makes sure that cells have the correct number of vertexes and that those vertexes are valid (on the master point list)
          - This check also makes sure cells have no duplicate vertexes.
          - In general this check is pretty cheap, and catches a great many common mesh problems.
          - This is the only check required for cell_type_t::POINT and cell_type_t::SEG cells
        - chk_cell_dimension:
          - This check makes sure that PYRAMIDS & HEXAHEDRONS are not co-planar and that TRIANGLES & QUADS are not co-linear.
          - Combined with chk_cell_vertexes this check is all that is required for cell_type_t::TRIANGLE cells
        - chk_cell_edges:
          - This check makes sure that no cell edges intersect in disallowed ways
          - This check when combined with the above checks, usually produce cells good enough for most software packages.
            - cell_type_t::QUAD cells might be concave or non-plainer.
            - cell_type_t::HEXAHEDRON cells might be concave or degenerate (bad edge-face intersections)
            - cell_type_t::PYRAMID cells might be concave.

      Levels Of Cell Goodness

       \verbatim                                                                                                              
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       |             | uniq_points       | uniq_points        | uniq_points        | uniq_points                   | uniq_points                   |
       |             | uniq_cells        | uniq_cells         | uniq_cells         | uniq_cells                    | uniq_cells                    |
       |             | chk_cell_vertexes | chk_cell_vertexes  | chk_cell_vertexes  | chk_cell_vertexes             | chk_cell_vertexes             |
       |             |                   | chk_cell_dimension | chk_cell_dimension | chk_cell_dimension            | chk_cell_dimension            |
       |             |                   |                    | chk_cell_edges     | chk_cell_edges                | chk_cell_edges                |
       |             |                   |                    |                    | check_cell_face_intersections | check_cell_face_intersections |
       |             |                   |                    |                    |                               | check_cell_faces_plainer      |
       |             |                   |                    |                    |                               | check_cell_convex             |
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       | vertex      | Perfect           |                    |                    |                               |                               |
       | segment     | Perfect           |                    |                    |                               |                               |
       | triangle    |                   | Perfect            |                    |                               |                               |
       | quad        |                   |                    | Good Enough        |                               | Perfect                       |
       | hexahedron  |                   |                    | Mostly OK          | Good Enough                   | Perfect                       |
       | pyramid     |                   |                    | Good Enough        |                               | Perfect                       |
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       \endverbatim


        @tparam uniq_points        Only allow unique points on the master point list
        @tparam uniq_cells         Only allow unique cells on the master cell list.  Won't prevent cells that are sub-cells of existing cells.
        @tparam chk_cell_vertexes  Do cell vertex checks (See: cell_stat_t)
        @tparam chk_cell_dimension Do cell dimension checks (See: cell_stat_t)
        @tparam chk_cell_edges     Do cell edge checks (See: cell_stat_t)
        @tparam eps                Epsilon used to detect zero */
  template <bool uniq_points, 
            bool uniq_cells, 
            bool chk_cell_vertexes,
            bool chk_cell_dimension,
            bool chk_cell_edges,
            double eps>
  class MR_cell_cplx {

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Master point list.

          cell complex objects have a master list of points.  Cells & data elements defined by referencing point indexes on this master list.  

          Points consist of three double values, think (x, y, z) in R^3.  Indexes in this list are used to identify points.  The first point added gets index
          0, and each successive point gets the next integer. */
      //@{
      /** Base field type for the vector space we use for points. */
      typedef double pnt_crd_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Hold a tuple of real values defining a point. */
      typedef std::array<pnt_crd_t, 3> pnt_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Integer type used to indentify/index points. */
      typedef int pnt_idx_t;
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Master point list. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** The index of the last point added via the add_point() method. */
      pnt_idx_t last_point_idx = -1;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** True if the last point given to the add_point() method was new -- i.e. not on the master point list. 
          Only updated if uniq_points is true.  See: last_point_added_was_new() */
      bool last_point_new = true;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Less operator for pnt_lst_uniq_t.
          @return If a & b are close in space, then return false.  Otherwise uses lexicographic ordering. 
          @bug Assumes 3==3. */
      struct pnt_less {
          bool operator()(const pnt_t& a, const pnt_t& b) const { return (((std::abs(a[0]-b[0]) > eps) ||
                                                                           (std::abs(a[1]-b[1]) > eps) ||
                                                                           (std::abs(a[2]-b[2]) > eps)) &&
                                                                          ((a[0] < b[0]) || 
                                                                           ((a[0] == b[0]) && 
                                                                            ((a[1] < b[1]) || 
                                                                             ((a[1] == b[1]) && (a[2] < b[2])))))); } 
      };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** List of points */
      typedef std::vector<pnt_t> pnt_idx_to_pnt_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Used to uniquify points and assign point index values */
      // typedef std::map<pnt_t, pnt_idx_t> pnt_to_pnt_idx_map_t;
      typedef std::map<pnt_t, pnt_idx_t, pnt_less> pnt_to_pnt_idx_map_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Maps points to point index -- used to detect physically identical points in R^3 */
      pnt_to_pnt_idx_map_t pnt_to_pnt_idx_map;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Maps point index to points -- the master point list */
      pnt_idx_to_pnt_t pnt_idx_to_pnt;
      //@}


    public:
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name 3D Vector Computations. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the magnitude */
      inline double vec3_two_norm(const pnt_t& pnt) const {
        return std::sqrt(vec3_self_dot_product(pnt));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the self dot product -- i.e. magnitude squared */
      inline double vec3_self_dot_product(const pnt_t& pnt) const {
        double tmp = 0.0;
        for(int i=0; i<3; ++i)
          tmp += pnt[i]*pnt[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the cross prodcut */
      inline double vec3_dot_product(const pnt_t& pnt1, const pnt_t& pnt2) const {
        double tmp = 0.0;
        for(int i=0; i<3; ++i)
          tmp += pnt1[i]*pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the cross prodcut */
      inline pnt_t vec3_cross_product(const pnt_t& pnt1, const pnt_t& pnt2) const {
        pnt_t tmp;
        tmp[0] = pnt1[1]*pnt2[2]-pnt1[2]*pnt2[1];
        tmp[1] = pnt1[2]*pnt2[0]-pnt1[0]*pnt2[2];
        tmp[2] = pnt1[0]*pnt2[1]-pnt1[1]*pnt2[0];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the vector diffrence */
      inline pnt_t vec3_diff(const pnt_t& pnt1, const pnt_t& pnt2) const {
        pnt_t tmp;
        for(int i=0; i<3; ++i)
          tmp[i] = pnt1[i] - pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the vector diffrence */
      inline double vec3_scalar_triple_product(const pnt_t& pnt1, const pnt_t& pnt2, const pnt_t& pnt3) const {
        return vec3_dot_product(pnt1, vec3_cross_product(pnt2, pnt3));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unitize the given point in place.  Return true if the result is valid, and false otherwise */
      inline bool vec3_unitize(pnt_t& pnt) const {
        double length = vec3_two_norm(pnt);
        if (std::abs(length) > eps) {
          for(int i=0; i<3; ++i) 
            pnt[i] = pnt[i]/length;
          return true;
        } else {
          return false;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the linear combination */
      inline pnt_t vec3_linear_combination(double scl1, const pnt_t& pnt1, double scl2, const pnt_t& pnt2) const {
        pnt_t tmp;
        for(int i=0; i<3; ++i)
          tmp[i] = scl1 * pnt1[i] + scl2 * pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Determinant of 3x3 matrix with given vectors as columns/rows. */
      inline double vec3_det3(const pnt_t& pnt1, const pnt_t& pnt2, const pnt_t& pnt3) const {
        //  MJR TODO NOTE <2024-07-22T15:48:31-0500> vec3_det3: UNTESTED! UNTESTED! UNTESTED! UNTESTED! 
        return (pnt1[0] * pnt2[1] * pnt3[2] - 
                pnt1[0] * pnt2[2] * pnt3[1] - 
                pnt1[1] * pnt2[0] * pnt3[2] + 
                pnt1[1] * pnt2[2] * pnt3[0] + 
                pnt1[2] * pnt2[0] * pnt3[1] - 
                pnt1[2] * pnt2[1] * pnt3[0]);
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name 3D Geometry Computations. 

       Methods with geomi_ prefix take point index types while methods with a geomr_ prefix take double tuples.  */
      //@{

      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** list of points */
      typedef std::vector<pnt_idx_t> pnt_idx_list_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Segment intersection types */
      enum class seg_isect_t { C0_EMPTY,        //!< Cardinality=0.        Intersection is the empty set.
                               C1_VERTEX1,      //!< Cardinality=1.        Intersection is a single, shared vertex.
                               C1_INTERIOR,     //!< Cardinality=1.        Intersection is a single point in interior of at least one segment.
                               CI_VERTEX2,      //!< Cardinality=infinite. Intersection is a segment -- equal to input segments.
                               CI_VERTEX1,      //!< Cardinality=infinite. Intersection is a segment -- exactially one vertex in the intersection.
                               CI_VERTEX0,      //!< Cardinality=infinite. Intersection is a segment -- no vertexes included
                               BAD_SEGMENT,     //!< At least one of the segments was degenerate
                             };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert an seg_isect_t to a string */
      std::string seg_isect_to_string(seg_isect_t seg_isect) const {
        switch(seg_isect) {
          case seg_isect_t::C0_EMPTY:     return std::string("C0_EMPTY");    break; 
          case seg_isect_t::C1_VERTEX1:   return std::string("C1_VERTEX1");  break; 
          case seg_isect_t::C1_INTERIOR:  return std::string("C1_INTERIOR"); break; 
          case seg_isect_t::CI_VERTEX2:   return std::string("CI_VERTEX2");  break; 
          case seg_isect_t::CI_VERTEX1:   return std::string("CI_VERTEX1");  break; 
          case seg_isect_t::CI_VERTEX0:   return std::string("CI_VERTEX0");  break; 
          case seg_isect_t::BAD_SEGMENT:  return std::string("BAD_SEGMENT"); break; 
        }
        return std::string(""); // Never get here, but some compilers can't figure that out. ;)
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Determine the nature of the intersection between two line segments */
      seg_isect_t geomi_seg_isect_type(pnt_idx_t ilin1pnt1, pnt_idx_t ilin1pnt2, pnt_idx_t ilin2pnt1, pnt_idx_t ilin2pnt2) const {
        //  MJR TODO NOTE geomi_seg_isect_type: Repeated point look-up slows things down
        //  MJR TODO NOTE geomi_seg_isect_type: Optimize
        //  MJR TODO NOTE geomi_seg_isect_type: Add unit tests for each code branch
        // Check for degenerate segments
        if (ilin1pnt1 == ilin1pnt2)
          return seg_isect_t::BAD_SEGMENT;
        if (ilin2pnt1 == ilin2pnt2)
          return seg_isect_t::BAD_SEGMENT;
        // Count unique points & break into cases
        std::set<pnt_idx_t> points_sorted;
        points_sorted.insert(ilin1pnt1);
        points_sorted.insert(ilin1pnt2);
        points_sorted.insert(ilin2pnt1);
        points_sorted.insert(ilin2pnt2);
        if (points_sorted.size() == 4) { // ...................................................... REMAINING CASES: C0_EMPTY, C1_INTERIOR, CI_VERTEX0
          if (geomi_pts_colinear(ilin1pnt1, ilin1pnt2, ilin2pnt1, ilin2pnt2)) { // ............... REMAINING CASES: C0_EMPTY,              CI_VERTEX0
            if ( (geomi_pnt_line_distance(ilin1pnt1, ilin1pnt2, ilin2pnt1, true) < eps) ||  
                 (geomi_pnt_line_distance(ilin1pnt1, ilin1pnt2, ilin2pnt2, true) < eps) ||
                 (geomi_pnt_line_distance(ilin2pnt1, ilin2pnt2, ilin1pnt1, true) < eps) ||  
                 (geomi_pnt_line_distance(ilin2pnt1, ilin2pnt2, ilin1pnt2, true) < eps) ) { // .. REMAINING CASES: CI_VERTEX0
              return seg_isect_t::CI_VERTEX0;
            } else { // ......................................................................... REMAINING CASES: C0_EMPTY
              return seg_isect_t::C0_EMPTY;
            }
          } else { // ........................................................................... REMAINING CASES: C0_EMPTY, C1_INTERIOR
            if (geomi_seg_isect1(ilin1pnt1, ilin1pnt2, ilin2pnt1, ilin2pnt2)) { // .............. REMAINING CASES: C1_INTERIOR
              return seg_isect_t::C1_INTERIOR;
            } else { // ......................................................................... REMAINING CASES: C0_EMPTY
              return seg_isect_t::C0_EMPTY;
            }
          }
          return seg_isect_t::C0_EMPTY;
        } else if (points_sorted.size() == 3) { // .............................................. REMAINING CASES: C1_VERTEX1, CI_VERTEX1          
          pnt_idx_t ipnt1, ipnt2, ipntc;
          if (ilin1pnt1 == ilin2pnt1) {
            ipntc = ilin1pnt1; ipnt1 = ilin1pnt2; ipnt2 = ilin2pnt2;
          } else if (ilin1pnt1 == ilin2pnt2) {
            ipntc = ilin1pnt1; ipnt1 = ilin1pnt2; ipnt2 = ilin2pnt1;
          } else if (ilin1pnt2 == ilin2pnt1) {
            ipntc = ilin1pnt2; ipnt1 = ilin1pnt1; ipnt2 = ilin2pnt2;
          } else if (ilin1pnt2 == ilin2pnt2) {
            ipntc = ilin1pnt2; ipnt1 = ilin1pnt1; ipnt2 = ilin2pnt1;
          } else { // Never get here.  Silences compiler warnings.
            ipntc = 0;         ipnt1 = 0;         ipnt2 = 0; 
          }
          if (geomi_pts_colinear(ipnt1, ipnt2, ipntc)) { // ..................................... REMAINING CASES: C1_VERTEX1, CI_VERTEX1          
            if ( (geomi_pnt_line_distance(ipnt1, ipntc, ipnt2, true) < eps) ||  
                 (geomi_pnt_line_distance(ipnt2, ipntc, ipnt1, true) < eps) ) { // .............. REMAINING CASES: CI_VERTEX1
              return seg_isect_t::CI_VERTEX1;
            } else { // ......................................................................... REMAINING CASES: CI_VERTEX1
              return seg_isect_t::C1_VERTEX1;
            }
          } else { // ........................................................................... REMAINING CASES: C1_VERTEX1
            return seg_isect_t::C1_VERTEX1;
          }
        } else if (points_sorted.size() == 2) { // .............................................. REMAINING CASES: CI_VERTEX2
          return seg_isect_t::CI_VERTEX2;
        } else { // (points_sorted.size() == 1) which is an error. // ........................... REMAINING CASES: CI_VERTEX2
          return seg_isect_t::CI_VERTEX2;
        }
        return seg_isect_t::C0_EMPTY;  // Should never get here...
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if two line segments intersect in a single point */
      bool geomi_seg_isect1(pnt_idx_t ilin1pnt1, pnt_idx_t ilin1pnt2, pnt_idx_t ilin2pnt1, pnt_idx_t ilin2pnt2) const {
        return geomr_seg_isect1(pnt_idx_to_pnt[ilin1pnt1], pnt_idx_to_pnt[ilin1pnt2], pnt_idx_to_pnt[ilin2pnt1], pnt_idx_to_pnt[ilin2pnt2]);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if two line segments intersect in a single point */
      bool geomr_seg_isect1(const pnt_t& lin1pnt1, const pnt_t& lin1pnt2, const pnt_t& lin2pnt1, const pnt_t& lin2pnt2) const {
        double denom = 
          lin1pnt1[0] * lin2pnt1[1] - lin1pnt1[0] * lin2pnt2[1] - lin1pnt1[1] * lin2pnt1[0] + lin1pnt1[1] * lin2pnt2[0] - 
          lin1pnt2[0] * lin2pnt1[1] + lin1pnt2[0] * lin2pnt2[1] + lin1pnt2[1] * lin2pnt1[0] - lin1pnt2[1] * lin2pnt2[0];
        if (std::abs(denom) < eps) // Lines are parallel
          return false;
        double numera = 
          lin1pnt1[0]*lin2pnt1[1] - lin1pnt1[0]*lin2pnt2[1] - 
          lin1pnt1[1]*lin2pnt1[0] + lin1pnt1[1]*lin2pnt2[0] + 
          lin2pnt1[0]*lin2pnt2[1] - lin2pnt1[1]*lin2pnt2[0];
        double numerb = 
          -(lin1pnt1[0]*lin1pnt2[1] - lin1pnt1[0]*lin2pnt1[1] - 
            lin1pnt1[1]*lin1pnt2[0] + lin1pnt1[1]*lin2pnt1[0] + 
            lin1pnt2[0]*lin2pnt1[1] - lin1pnt2[1]*lin2pnt1[0]);
        double ua = numera/denom;
        double ub = numerb/denom;
        if ( (ua < 0) || (ub < 0) || (ua > 1) || (ub > 1) )
          return false;
        double eq3 = numera/denom * (lin1pnt2[2] - lin1pnt1[2]) - numerb/denom * (lin2pnt2[2] - lin2pnt1[2]) + lin2pnt1[2] + lin1pnt1[2];
        if (std::abs(eq3) < eps) // Equation in third coordinate is satisfied
          return true;
        else
          return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between a point and a line. 
          See: geomr_pnt_line_distance(). */
      double geomi_pnt_line_distance(pnt_idx_t ilinpnt1, pnt_idx_t ilinpnt2, pnt_idx_t ipnt, bool seg_distance) const {
        return geomr_pnt_line_distance(pnt_idx_to_pnt[ilinpnt1], pnt_idx_to_pnt[ilinpnt2], pnt_idx_to_pnt[ipnt], seg_distance);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between a point and a line.
          The result depends on the value of seg_distance:
           - seg_distance==false: Distance between the point pnt and the line containing linpnt1 & linpnt2.
           - seg_distance==true:  Distance between the point pnt and the line segment with endpoints linpnt1 & linpnt2. */
      double geomr_pnt_line_distance(const pnt_t& linpnt1, const pnt_t& linpnt2, const pnt_t& pnt, bool seg_distance) const {
        double segd = geomr_pnt_pnt_distance(linpnt1, linpnt2);
        pnt_t d, p;
        double t = 0;
        for(int i=0; i<3; i++) {
          d[i] = (linpnt2[i] - linpnt1[i]) / segd;
          t += (pnt[i] - linpnt1[i])*d[i];
        }
        for(int i=0; i<3; i++) {
          p[i] = linpnt1[i] + t * d[i];
        }
        // p is the point on the line nearest pnt
        if (seg_distance) {
          double dp1 = geomr_pnt_pnt_distance(linpnt1, p);
          double dp2 = geomr_pnt_pnt_distance(linpnt2, p);
          if (std::abs((dp1+dp2)-segd) > eps) // MJR TODO NOTE: check logic -- dp1>segd || dp2>segd?
            return std::min(geomr_pnt_pnt_distance(linpnt1, pnt), geomr_pnt_pnt_distance(linpnt2, pnt));
        }
        return geomr_pnt_pnt_distance(p, pnt);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the euculedian (2 norm) distance between two points */
      double geomr_pnt_pnt_distance(const pnt_t& pnt1, const pnt_t& pnt2) const {
        return vec3_two_norm(vec3_diff(pnt1, pnt2));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the normal of a triangle (or a plane defined by 3 points) */
      inline pnt_t geomr_tri_normal(const pnt_t& tripnt1, const pnt_t& tripnt2, const pnt_t& tripnt3, bool unit) const {
        pnt_t basisv1 = vec3_diff(tripnt1, tripnt2);  // basis vectors for pln containing triagnel
        pnt_t basisv2 = vec3_diff(tripnt3, tripnt2);  // basis vectors for pln containing triagnel
        pnt_t normal = vec3_cross_product(basisv1, basisv2); // normal vector for tri. n=pld1xpld2
        if (unit)
          vec3_unitize(normal);
        return normal;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the distance between a point and a plane */
      double geomr_pnt_pln_distance(const pnt_t& plnpnt1, const pnt_t& plnpnt2, const pnt_t& plnpnt3, const pnt_t& pnt) const {
        pnt_t n = geomr_tri_normal(plnpnt1, plnpnt2, plnpnt3, true);
        return std::abs(vec3_dot_product(n, pnt) - vec3_dot_product(n, plnpnt2));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the distance between a point and a triangle */
      double geomr_pnt_tri_distance(const pnt_t& tripnt1, const pnt_t& tripnt2, const pnt_t& tripnt3, const pnt_t& pnt) const {
        //  MJR TODO NOTE <2024-07-22T15:48:31-0500> geomr_pnt_tri_distance: UNTESTED! UNTESTED! UNTESTED! UNTESTED! 
        pnt_t basisv1 = vec3_diff(tripnt1, tripnt2);  // basis vectors for pln containing triagnel
        pnt_t basisv2 = vec3_diff(tripnt3, tripnt2);  // basis vectors for pln containing triagnel
        pnt_t normal = vec3_cross_product(basisv1, basisv2); // normal vector for tri. ax+by+cz+d=0, a=normal[0], b=normal[1], c=normal[2]
        vec3_unitize(normal); 
        double d = -vec3_dot_product(normal, tripnt2);            // ax+by+cz+d=0
        double lambda = vec3_dot_product(normal, pnt) + d;
        pnt_t q = vec3_linear_combination(1.0, pnt, lambda, normal); // q is the point in the plane closest to pnt
        double denom =  basisv1[0] * basisv2[1] - basisv2[0] * basisv1[1]; // If zero, then triangle is broken!
        double u     = (q[0]       * basisv2[1] - basisv2[0] *       q[1]) / denom;
        double v     = (basisv1[0] *       q[1] -       q[0] * basisv1[1]) / denom;
        double dd    = std::abs(u*basisv1[2] + v*basisv2[2] - q[2]);
        if ( (u>=0) && (v>=0) && ((u+v)<=1) && (dd<eps) ) { // q is in the triangle =>  Use the plane distance
          return std::abs(lambda);
        } else {                                            // q is not in the triangle =>  Distance will be minimum distance to an edge.
          double d1 = geomr_pnt_line_distance(tripnt1, tripnt2, pnt, true);
          double d2 = geomr_pnt_line_distance(tripnt2, tripnt3, pnt, true);
          double d3 = geomr_pnt_line_distance(tripnt3, tripnt1, pnt, true);
          return std::min(std::min(d1, d2), d3);
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points (just one point for this function) are zero.
          Checks if the infinity norm is less than epsilon.*/
      bool geomr_pnt_zero(const pnt_t& p1) const {
        return ((std::abs(p1[0]) < eps) &&
                (std::abs(p1[1]) < eps) &&
                (std::abs(p1[2]) < eps));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear */
      bool geomi_pts_colinear(pnt_idx_t pi1, pnt_idx_t pi2, pnt_idx_t pi3, pnt_idx_t pi4) const {
        return ( geomr_pts_colinear(pnt_idx_to_pnt[pi1], pnt_idx_to_pnt[pi2], pnt_idx_to_pnt[pi3]) &&
                 geomr_pts_colinear(pnt_idx_to_pnt[pi1], pnt_idx_to_pnt[pi2], pnt_idx_to_pnt[pi4]) );
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear. */
      bool geomi_pts_colinear(pnt_idx_t pi1, pnt_idx_t pi2, pnt_idx_t pi3) const {
        return geomr_pts_colinear(pnt_idx_to_pnt[pi1], pnt_idx_to_pnt[pi2], pnt_idx_to_pnt[pi3]);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear */
      bool geomr_pts_colinear(pnt_t p1, pnt_t p2, pnt_t p3) const {
        return geomr_pnt_zero(vec3_cross_product(vec3_diff(p1, p2), vec3_diff(p1, p3)));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      bool geomi_pts_coplanar(const pnt_idx_list_t& pnt_list) const {
        if (pnt_list.size() > 3) {
          if ( !(geomi_pts_coplanar(pnt_list[0], pnt_list[1], pnt_list[2], pnt_list[3])))
            return false;
          for(decltype(pnt_list.size()) i=4; i<pnt_list.size(); i++)
            if ( !(geomi_pts_coplanar(pnt_list[0], pnt_list[1], pnt_list[2], pnt_list[i])))
              return false;
        }
        return true;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      bool geomi_pts_coplanar(pnt_idx_t pi1, pnt_idx_t pi2, pnt_idx_t pi3, pnt_idx_t pi4) const {
        return geomr_pts_coplanar(pnt_idx_to_pnt[pi1], pnt_idx_to_pnt[pi2], pnt_idx_to_pnt[pi3], pnt_idx_to_pnt[pi4]);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      bool geomr_pts_coplanar(pnt_t p1, pnt_t p2, pnt_t p3, pnt_t p4) const {
        return (std::abs(vec3_scalar_triple_product(vec3_diff(p3, p1), vec3_diff(p2, p1), vec3_diff(p4, p1))) < eps);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute intersection of a segment and a triangle */
      bool geomr_seg_tri_intersection(pnt_t tp1, pnt_t tp2, pnt_t tp3, pnt_t sp1, pnt_t sp2) const {
        //  MJR TODO NOTE <2024-07-11T16:08:16-0500> geomr_seg_tri_intersection: implement
        //  MJR TODO NOTE <2024-07-11T16:08:27-0500> geomr_seg_tri_intersection: Should this be a bool or an enum?
        return (tp1[0]+tp2[0]+tp3[0]+sp1[0]+sp2[0]>1.0);
      }
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Class utilizes. */
      //@{
      void clear() {
        last_point_idx = -1;
        last_point_new = true;
        pnt_to_pnt_idx_map.clear();
        pnt_idx_to_pnt.clear();
        pdata_sdat.clear();
        pdata_vdat.clear();
        cell_lst.clear();
        uniq_cell_lst.clear();
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Master point list. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns the index of the last point given to the add_point() method. */
      pnt_idx_t idx_of_last_point_added() const {
        return last_point_idx;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns true if the last point given to the add_point() method was a new point. */
      pnt_idx_t last_point_added_was_new() const {
        return last_point_new;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a point.
          Cases
           - Any of the coordinate values in the given point are NaN: last_point_idx=-1, last_point_new=false.
           - The given point is already on the list: last_point_idx is ste to the existing point's index, and last_point_new=false
           - The given point is not on the list: last_point_idx is set to th enew piont's index, and last_point_new=true
          Note that last_point_idx is always the resturn value. */
      pnt_idx_t add_point(pnt_t new_pnt) {
        if (std::isnan(new_pnt[0]) || std::isnan(new_pnt[0]) || std::isnan(new_pnt[0])) {
          last_point_idx = -1;
          last_point_new = false;
        } else {
          if constexpr (uniq_points) {
            if (pnt_to_pnt_idx_map.contains(new_pnt)) {
              /* Point is already in list */
              last_point_idx = pnt_to_pnt_idx_map[new_pnt];
              last_point_new = false;
            } else {
              /* Point is not already in list */
              last_point_idx = static_cast<pnt_idx_t>(pnt_idx_to_pnt.size());
              pnt_to_pnt_idx_map[new_pnt] = last_point_idx;
              pnt_idx_to_pnt.push_back(new_pnt);
              last_point_new = true;
            }
          } else {
            last_point_idx = static_cast<pnt_idx_t>(pnt_idx_to_pnt.size()); 
            pnt_idx_to_pnt.push_back(new_pnt);
            last_point_new = true;
          }
        }
        return last_point_idx;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Print number of points in master point list.
          Note the return type is pnt_idx_t (a signed integer type) and not a size_t. */
      pnt_idx_t num_points() const {
        return static_cast<pnt_idx_t>(pnt_idx_to_pnt.size());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a pnt_t to a string representation */
      std::string pnt_to_string(pnt_t x) const {
        std::ostringstream convert;
        convert << "[ ";
        for(auto c: x)
          convert << std::setprecision(5) << static_cast<double>(c) << " ";
        convert << "]";
        return(convert.str());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Print all points to STDOUT. 

          @param max_num_print Maximum number of points to print.  Use 0 to print all points. */
      void print_all_points(int max_num_print) const {
        int numPrinted = 0;
        if (num_points() > 0) {
          std::cout << "POINTS BEGIN (" << num_points() << ")" << std::endl;
          for(pnt_idx_t pnt_idx = 0; pnt_idx<num_points(); ++pnt_idx) {
            std::cout << "  " << pnt_idx << ": " << pnt_to_string(pnt_idx_to_pnt[pnt_idx]) << std::endl;
            numPrinted++;
            if ((max_num_print > 0) && (numPrinted >= max_num_print)) {
              std::cout << "  Maximum number of points reached.  Halting tree dump." << std::endl;
              break;
            }
          }
          std::cout << "POINTS END" << std::endl;
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Point data sets.

          Cell complex objects may have data sets associated with points.  Each dataset has a name (an std::string), and a type (scalar or 3D vector).  */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      typedef std::string pdata_name_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a single scalar data value for a single point */
      typedef pnt_crd_t sdat_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a single vector data value for a single point */
      typedef pnt_t vdat_t;
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Point data sets. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a scalar value for each point in the master point list. */
      typedef std::vector<pnt_crd_t> pnt_idx_to_sdat_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a group of named scalar data sets. */
      typedef std::map<pdata_name_t, pnt_idx_to_sdat_t> dat_name_to_sdat_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a 3D vector value for each point in the master point list. */
      typedef std::vector<pnt_t> pnt_idx_to_vdat_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a group of named vector data sets. */
      typedef std::map<pdata_name_t, pnt_idx_to_vdat_t> dat_name_to_vdat_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** House all scalar data sets. */
      dat_name_to_sdat_lst_t pdata_sdat;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** House all vector data sets. */
      dat_name_to_vdat_lst_t pdata_vdat;
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Point data sets. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a scalar point data value */
      void add_data(pdata_name_t pdata_name, pnt_idx_t pnt_idx, sdat_t da_data) {
        int num_dat = static_cast<int>(pdata_sdat[pdata_name].size());
        if (pnt_idx >= num_dat)
          pdata_sdat[pdata_name].resize(pnt_idx+1);
        pdata_sdat[pdata_name][pnt_idx] = da_data;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a vector point data value */
      void add_data(pdata_name_t pdata_name, pnt_idx_t pnt_idx, vdat_t da_data) {
        int num_dat = static_cast<int>(pdata_vdat[pdata_name].size());
        if (pnt_idx >= num_dat)
          pdata_vdat[pdata_name].resize(pnt_idx+1);
        pdata_vdat[pdata_name][pnt_idx] = da_data;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a vector point data value to last point added if the point was new */
      void add_data_if_new(pdata_name_t pdata_name, vdat_t da_data) {
        if (last_point_new)
          add_data(pdata_name, last_point_idx, da_data);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a scalar point data value to last point added if the point was new */
      void add_data_if_new(pdata_name_t pdata_name, sdat_t da_data) {
        if (last_point_new)
          add_data(pdata_name, last_point_idx, da_data);
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell Status.  */
      enum class cell_stat_t { GOOD,           //!< Looks like a good cell                                       N/A
                               TOO_FEW_PNT,    //!< List of points was empty                                     check_cell_vertexes
                               TOO_MANY_PNT,   //!< List of points was too long (>5)                             check_cell_vertexes
                               NEG_PNT_IDX,    //!< Negative point index (i.e. not on the master point list)     check_cell_vertexes
                               BIG_PNT_IDX,    //!< Point index was too big (i.e. not on the master point list)  check_cell_vertexes
                               DUP_PNT,        //!< At least two points have the same index                      check_cell_vertexes
                               DIM_LOW,        //!< Dimension low (degenerate)                                   check_cell_dimension
                               BAD_EDGEI,      //!< Bad edge-edge intersection                                   check_cell_edge_intersections
                               BAD_FACEI,      //!< Bad face-edge intersection                                   check_cell_face_intersections
                               FACE_BENT,      //!< A face (QUAD, HEXAHEDRON, or PYRAMID) was not plainer        check_cell_faces_plainer
                               CONCAVE,        //!< (QUAD, HEXAHEDRON, or PYRAMID) was concave                   TBD
                             };

      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell Types. */
      enum class cell_type_t { POINT,
                               SEGMENT,
                               TRIANGLE,
                               QUAD,
                               HEXAHEDRON,
                               PYRAMID,
                             };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a poly cell -- a list of point indexes */
      typedef pnt_idx_list_t cell_t;
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type for a list of poly cells */
      typedef std::vector<cell_t> cell_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** House all poly cells */
      cell_lst_t cell_lst;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Used to uniquify cells */
      // typedef std::map<pnt_t, pnt_idx_t> pnt_to_pnt_idx_map_t;
      typedef std::set<cell_t> uniq_cell_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unique cell list. */
      uniq_cell_lst_t uniq_cell_lst;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** True if the last cell given to the add_cell() method was new -- i.e. not on the master cell list. 
          Only updated if uniq_cells is true.  Value is invalid if last_cell_stat is NOT cell_stat_t::GOOD. 
          See: last_cell_added_was_new() */
      bool last_cell_new = true;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Status of the last cell added via add_cell() method. 
       Only updated if (chk_cell_vertexes || chk_cell_dimension | chk_cell_edges) is true. See: status_of_last_cell_added() */
      cell_stat_t last_cell_stat = cell_stat_t::GOOD;
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell segment/face/etc structure.  */
      typedef cell_lst_t cell_structure_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Number of cells. */
      inline int num_cells() const {
        return static_cast<int>(cell_lst.size());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_type_t value to a list of vertexes. */
      inline cell_structure_t& cell_type_to_structure(cell_type_t cell_type, int dimension) const {
        //  MJR TODO NOTE <2024-07-11T16:06:42-0500> cell_type_to_structure: make sure polygons are all oriented correctly
        int logical_dim = cell_type_to_dimension(cell_type);
        if ( (dimension < 0) || (dimension > logical_dim) )
          dimension = logical_dim;
        if (logical_dim > 3) {
          std::cerr << "MR_cell_cplx API USAGE ERROR: Maximum supported dimension is 3" << std::endl;
          exit(0);
        }
        int idx = 0;
        switch(cell_type) {
          case cell_type_t::POINT:       idx = 0; break;
          case cell_type_t::SEGMENT:     idx = 1; break;
          case cell_type_t::TRIANGLE:    idx = 2; break;
          case cell_type_t::QUAD:        idx = 3; break;
          case cell_type_t::PYRAMID:     idx = 4; break;
          case cell_type_t::HEXAHEDRON:  idx = 5; break;
        }
        static std::vector<std::vector<cell_structure_t>> cst = {{{{0}},                             // vertex
                                                                  {{0},{1}},                         // segment
                                                                  {{0},{1},{2}},                     // triangle
                                                                  {{0},{1},{2},{3}},                 // Quad
                                                                  {{0},{1},{2},{3},                  // Pyramid: Base vertexes
                                                                   {4}},                             // Pyramid: Tip vertex
                                                                  {{0},{1},{2},{3},                  // Hexahedron: Back vertexes
                                                                   {4},{5},{6},{7}}},                // Hexahedron: Front vertexes
                                                                 {{},                                // vertex      
                                                                  {{0,1}},                           // segment     
                                                                  {{0,1},{1,2},{2,0}},               // triangle    
                                                                  {{0,1},{1,2},{2,3},{3,0}},         // Quad        
                                                                  {{0,1},{1,2},{2,3},{3,0},          // Pyramid: Base segments
                                                                   {0,4},{1,4},{2,4},{3,4}},         // Pyramid: Side sets
                                                                  {{0,1},{1,2},{2,3},{3,0},          // Hexahedron: Back segments
                                                                   {4,5},{5,6},{6,7},{7,4},          // Hexahedron: Front segments
                                                                   {0,4},{1,5},{2,6},{3,7}}},        // Hexahedron: Back to front segments
                                                                 {{},                                // vertex      
                                                                  {},                                // segment     
                                                                  {{0,1,2}},                         // triangle    
                                                                  {{0,1,2,3}},                       // Quad        
                                                                  {{0,1,2,3},                        // Pyramid: Base face
                                                                   {0,1,4},                          // Pyramid: Back face
                                                                   {1,2,4},                          // Pyramid: Left face
                                                                   {2,3,4},                          // Pyramid: Front face
                                                                   {3,0,4}},                         // Pyramid: Right faces
                                                                  {{0,1,2,3},                        // Hexahedron: Back face
                                                                   {4,5,6,7},                        // Hexahedron: Front face
                                                                   {0,3,7,4},                        // Hexahedron: Left face
                                                                   {2,3,7,6},                        // Hexahedron: Top face
                                                                   {1,2,6,5},                        // Hexahedron: Right face
                                                                   {0,1,4,5}}},                      // Hexahedron: Bottom face
                                                                 {{},                                // vertex      
                                                                  {},                                // segment     
                                                                  {},                                // triangle    
                                                                  {},                                // Quad        
                                                                  {{0,1,2,3,4}},                     // Pyramid
                                                                  {{0,1,2,3,4,5,6,7}}}};             // Hexahedron
        return (cst[dimension][idx]);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_type_t value to the logical dimension of the cell.

          Note that the vertexes of a cell_type_t::QUAD cell might not be coplainer; however, the cell is logically a 2D cell.  Similarly, the vertexes of a
          cell_type_t::PYRAMID might be coplainer; however, the cell is logically a 3D cell.  */
      inline int cell_type_to_dimension(cell_type_t cell_type) const {
        switch(cell_type) {
          case cell_type_t::POINT:       return (0); break;
          case cell_type_t::SEGMENT:     return (1); break;
          case cell_type_t::TRIANGLE:    return (2); break;
          case cell_type_t::QUAD:        return (2); break;
          case cell_type_t::PYRAMID:     return (3); break;
          case cell_type_t::HEXAHEDRON:  return (3); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_type_t value to the number of points required for cell type. */
      inline int cell_type_to_req_pt_cnt(cell_type_t cell_type) const {
        switch(cell_type) {
          case cell_type_t::POINT:       return (1); break;
          case cell_type_t::SEGMENT:     return (2); break;
          case cell_type_t::TRIANGLE:    return (3); break;
          case cell_type_t::QUAD:        return (4); break;
          case cell_type_t::PYRAMID:     return (5); break;
          case cell_type_t::HEXAHEDRON:  return (8); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_type_t value to the VTK cell type integer. */
      inline int cell_type_to_vtk_type(cell_type_t cell_type) const {
        switch(cell_type) {
          case cell_type_t::POINT:       return ( 1); break;
          case cell_type_t::SEGMENT:     return ( 3); break;
          case cell_type_t::TRIANGLE:    return ( 5); break;
          case cell_type_t::QUAD:        return ( 9); break;
          case cell_type_t::HEXAHEDRON:  return (12); break;
          case cell_type_t::PYRAMID:     return (14); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_type_t value to the VTK cell type integer. */
      inline std::string cell_type_to_string(cell_type_t cell_type) const {
        switch(cell_type) {
          case cell_type_t::POINT:       return ("POINT"     ); break;
          case cell_type_t::SEGMENT:     return ("SEGMENT"   ); break;
          case cell_type_t::TRIANGLE:    return ("TRIANGLE"  ); break;
          case cell_type_t::QUAD:        return ("QUAD"      ); break;
          case cell_type_t::HEXAHEDRON:  return ("HEXAHEDRON"); break;
          case cell_type_t::PYRAMID:     return ("PYRAMID"   ); break;
        }
        return ""; // Never get here.  
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert number of points in a cell to cell_type_t. */
      inline cell_type_t req_pt_cnt_to_cell_type(std::vector<int>::size_type num_points) const {
        switch(num_points) {
          case 1: return (cell_type_t::POINT);      break;
          case 2: return (cell_type_t::SEGMENT);    break;
          case 3: return (cell_type_t::TRIANGLE);   break;
          case 4: return (cell_type_t::QUAD);       break;
          case 5: return (cell_type_t::PYRAMID);    break;
          case 8: return (cell_type_t::HEXAHEDRON); break;
        }
        return (cell_type_t::POINT); // Never get here
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t to a bool (true if GOOD). */
      inline bool cell_stat_is_good(cell_stat_t cell_stat) const {
        return (cell_stat == cell_stat_t::GOOD);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t to a bool (true if BAD). */
      inline bool cell_stat_is_bad(cell_stat_t cell_stat) const {
        return (cell_stat != cell_stat_t::GOOD);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t enum value to a string. */
      std::string cell_stat_to_string(cell_stat_t cell_stat) const {
        switch(cell_stat) {
          case cell_stat_t::GOOD:            return (std::string("GOOD"));            break;
          case cell_stat_t::TOO_FEW_PNT:     return (std::string("TOO_FEW_PNT"));     break;
          case cell_stat_t::TOO_MANY_PNT:    return (std::string("TOO_MANY_PNT"));    break;
          case cell_stat_t::NEG_PNT_IDX:     return (std::string("NEG_PNT_IDX"));     break;
          case cell_stat_t::BIG_PNT_IDX:     return (std::string("BIG_PNT_IDX"));     break;
          case cell_stat_t::DUP_PNT:         return (std::string("DUP_PNT"));         break;
          case cell_stat_t::DIM_LOW:         return (std::string("DIM_LOW"));         break;
          case cell_stat_t::BAD_EDGEI:       return (std::string("BAD_EDGEI"));       break;
          case cell_stat_t::BAD_FACEI:       return (std::string("BAD_FACEI"));       break;
          case cell_stat_t::FACE_BENT:       return (std::string("FACE_BENT"));       break;
          case cell_stat_t::CONCAVE:         return (std::string("CONCAVE"));         break;
        }
        return std::string("ERROR");
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Perform cell vertex checks

          @param cell_type The type for a_cell.
          @param a_cell    The cell to test.  */
      cell_stat_t check_cell_vertexes(cell_type_t cell_type, cell_t a_cell) const {
        // Check number of points
        std::vector<int>::size_type a_cell_len  = a_cell.size();
        std::vector<int>::size_type req_num_pts = cell_type_to_req_pt_cnt(cell_type);
        if (a_cell_len < req_num_pts)
          return cell_stat_t::TOO_FEW_PNT;
        if (a_cell_len > req_num_pts)
          return cell_stat_t::TOO_MANY_PNT;
        // Check for negative point index
        if (std::any_of(a_cell.cbegin(), a_cell.cend(), [](pnt_idx_t i) { return (i < 0); }))
          return cell_stat_t::NEG_PNT_IDX;
        // Check for too big point index
        if (std::any_of(a_cell.cbegin(), a_cell.cend(), [this](pnt_idx_t i) { return (i >= num_points()); }))
          return cell_stat_t::BIG_PNT_IDX;
        // Check for duplicate points
        if (a_cell_len > 1) {
          std::set<pnt_idx_t> a_cell_pnt_sorted;
          for(pnt_idx_t pnt_idx: a_cell) {
            if (a_cell_pnt_sorted.contains(pnt_idx))
              return cell_stat_t::DUP_PNT;          
            a_cell_pnt_sorted.insert(pnt_idx);
          }
        }
        // Return GOOD
        return cell_stat_t::GOOD;          
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Perform cell dimension (2D cells must not be colinear & 3D cells must not be coplainer) 

          @warning This function assumes that check_cell_vertexes() would have returned GOOD for the cell being checked!

          @warning Will not detect degenerate cell_type_t::SEGMENT as that is handled by check_cell_vertexes.

          @param cell_type The type for a_cell.
          @param a_cell    The cell to test.  */
      cell_stat_t check_cell_dimension(cell_type_t cell_type, cell_t a_cell) const {
        if (cell_type == cell_type_t::TRIANGLE) {
          if (geomi_pts_colinear(a_cell[0], a_cell[1], a_cell[2]))
            return cell_stat_t::DIM_LOW;
        } else if (cell_type == cell_type_t::QUAD) {
          if (geomi_pts_colinear(a_cell[0], a_cell[1], a_cell[2], a_cell[3]))
            return cell_stat_t::DIM_LOW;
        } else if (cell_type == cell_type_t::HEXAHEDRON) {
          if ( geomi_pts_coplanar(a_cell))
            return cell_stat_t::DIM_LOW;
        } else if (cell_type == cell_type_t::PYRAMID) {
          if ( geomi_pts_coplanar(a_cell))
            return cell_stat_t::DIM_LOW;
        }
        return cell_stat_t::GOOD;          
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Checks that cell edges have expected intersections.

          Checks that every pair of cell edges has the correct intersection type (i.e. intersect in a single vertex or an empty intersection).  If a bad
          intersection is detected, then a cell_stat_t of cell_stat_t::BAD_EDGEI will be returned.  Otherwise cell_stat_t::GOOD will be returned.  

          @warning This function assumes that check_cell_vertexes() would have returned GOOD for the cell being checked!

          Note that this function will detect some conditions caught by other checks.  For example, this function will return cell_stat_t::BAD_EDGEI for
          cell_type_t::SEG cells for which check_cell_dimension() returns cell_stat_t::DUP_PNT.  Note the cell_stat_t values are different depending upon which
          check function is called.
          
          @param cell_type The type for a_cell.
          @param a_cell    The cell to test.  */
      cell_stat_t check_cell_edge_intersections(cell_type_t cell_type, cell_t a_cell) const {
        cell_structure_t& segs = cell_type_to_structure(cell_type, 1);
        if ( !(segs.empty())) {
          for(decltype(segs.size()) i=0; i<segs.size()-1; i++) {
            for(decltype(segs.size()) j=i+1; j<segs.size(); j++) {
              //std::cout << segs[i][0] << "--" << segs[i][1] << " CAP " << segs[j][0] << "--" << segs[j][1] << std::endl;
              std::set<pnt_idx_t> points_sorted;
              points_sorted.insert(segs[i][0]);
              points_sorted.insert(segs[i][1]);
              points_sorted.insert(segs[j][0]);
              points_sorted.insert(segs[j][1]);
              auto it = geomi_seg_isect_type(a_cell[segs[i][0]], a_cell[segs[i][1]], a_cell[segs[j][0]], a_cell[segs[j][1]]);
              if(points_sorted.size() == 4) {
                if (it != seg_isect_t::C0_EMPTY)
                  return cell_stat_t::BAD_EDGEI;
              } else if(points_sorted.size() == 3) {
                if (it != seg_isect_t::C1_VERTEX1) 
                  return cell_stat_t::BAD_EDGEI;
              } else {
                return cell_stat_t::BAD_EDGEI;
              }
            }
          }
        }
        return cell_stat_t::GOOD;          
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Checks that cell faces have expected intersections.
          
          @param cell_type The type for a_cell.
          @param a_cell    The cell to test.  */
      cell_stat_t check_cell_face_intersections(cell_type_t cell_type, cell_t a_cell) const {
        //  MJR TODO NOTE check_cell_face_intersections: Implement
        if (cell_type == cell_type_t::HEXAHEDRON) {
          if ( geomi_pts_coplanar(a_cell))
            return cell_stat_t::DIM_LOW;
        } else if (cell_type == cell_type_t::PYRAMID) {
          if ( geomi_pts_coplanar(a_cell))
            return cell_stat_t::DIM_LOW;
        }
        return cell_stat_t::GOOD;          
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if the vertexes of a cell are coplainer.

          Returns cell_stat_t::FACE_BENT if the cell is *not* plainer, and cell_stat_t::GOOD if it is.

          Note that most visualization packages are tolerant of non-plainer faces, and can render them just fine.  For example, most applications will
          automatically split a non-plainer cell_type_t::QUAD into two triangles.  That said such faces can be an issue when using a tessellation for
          computation.  Such issues are most severe for things like FEM, but they can also show up for common visualization computations like the extraction
          of level curves/surfaces.

          @param cell_type The type for a_cell.
          @param a_cell    The cell to test.  */
      cell_stat_t check_cell_faces_plainer(cell_type_t cell_type, cell_t a_cell) const {
        const cell_structure_t& face_structures = cell_type_to_structure(cell_type, 2);
        for(auto face_structure: face_structures) {
          cell_t face;
          for(auto idx: face_structure) 
            face.push_back(a_cell[idx]);
          // std::cout << "[ ";
          // for(auto v: face) 
          //   std::cout << v << " ";
          // std::cout << "]" << std::endl;
          if ( !(geomi_pts_coplanar(face)))
            return cell_stat_t::FACE_BENT;
        }
        return cell_stat_t::GOOD;             
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add parts of a cell of the specified dimension.

          Example 1: a valid cell_type_t::PYRAMID cell added wth dimension == 2 results in three cell_type_t::TRIANGLE cells and one cell_type_t::QUAD cell.

          Example 2: a valid cell_type_t::PYRAMID cell added wth dimension == 3 results in one cell_type_t::PYRAMID cell added

          Example 3: a valid cell_type_t::PYRAMID cell added wth dimension == 1 results in 8 cell_type_t::SEG cells added

          @param cell_type The type of cell to add.
          @param new_cell  The cell to add. 
          @param dimension The dimension of the parts to add.
          @return Number of cells added */
      int add_cell(cell_type_t cell_type, cell_t new_cell, int dimension) {
        int num_added = 0;
        if ( (dimension < 0) || (dimension >= cell_type_to_dimension(cell_type)) ) {
          if (add_cell(cell_type, new_cell))
            num_added++;
        } else { // We need to break the cell up into lower dimensional bits, and add the bits.
          const cell_structure_t& cell_parts = cell_type_to_structure(cell_type, dimension);
          for(auto cell_part: cell_parts) {
            cell_t newer_cell;
            for(auto idx: cell_part) 
              newer_cell.push_back(new_cell[idx]);
            if (add_cell(req_pt_cnt_to_cell_type(newer_cell.size()), newer_cell))
              num_added++;
          }
        }
        return num_added;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a cell
          @param cell_type The type of cell to add.
          @param new_cell  The cell to add. 
          @return A boolean indicateing success
          @retval true  The cell was added or had been added previously
          @retval false The cell could not be added (because of a failed geometric check) */
      bool add_cell(cell_type_t cell_type, cell_t new_cell) {
        // Check vertexes if required
        if constexpr (chk_cell_vertexes) {
          last_cell_stat = check_cell_vertexes(cell_type, new_cell);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Check dimension if required
        if constexpr (chk_cell_dimension) {
          last_cell_stat = check_cell_dimension(cell_type, new_cell);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Check edges
        if constexpr (chk_cell_edges) {
          last_cell_stat = check_cell_edge_intersections(cell_type, new_cell);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Geom was good or we didn't need to check.  
        if constexpr (uniq_cells) {
          cell_t sorted_cell = new_cell;
          std::sort(sorted_cell.begin(), sorted_cell.end());
          if (uniq_cell_lst.contains(sorted_cell)) {
            last_cell_new = false;
          } else {
            last_cell_new = true;
            cell_lst.push_back(new_cell);
            uniq_cell_lst.insert(sorted_cell);
          }
        } else {
          cell_lst.push_back(new_cell);
        }
        return true;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Print all cells to STDOUT. 
          @param max_num_print Maximum number of cells to print.  Use 0 to print all cells. */
      void print_all_cells(int max_num_print) const {
        int numPrinted = 0;
        if (num_cells() > 0) {
          std::cout << "CELLS BEGIN (" << num_cells() << ")" << std::endl;
          for(int i=0; i<num_cells(); i++) {
            std::cout << "  ";
            for(auto& vert: cell_lst[i]) {
              std::cout << vert << " ";
            }
            std::cout << "   " << cell_type_to_string(req_pt_cnt_to_cell_type(cell_lst[i].size())) << std::endl;
            numPrinted++;
            if ((max_num_print > 0) && (numPrinted >= max_num_print)) {
              std::cout << "  Maximum number of cells reached.  Halting tree dump." << std::endl;
              break;
            }
          }
          std::cout << "CELLS END" << std::endl;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns the status of the last cell given to the add_cell() method.
          If (chk_cell_vertexes || chk_cell_dimension | chk_cell_edges) is true, this value is updated each time add_cell() is called.  
          Otherwise its value is always cell_stat_t::GOOD. */
      cell_stat_t status_of_last_cell_added() const {
        return last_cell_stat;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns true if the last cell given to the add_cell() was new and the return value of status_of_last_cell_added() is cell_stat_t::GOOD.
          If uniq_cells is true, this value is updated each time add_cell() is called.  Otherwise its value is always true. */
      bool last_cell_added_was_new() const {
        return last_cell_new;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name I/O. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to an XML VTK unstructured grid file

          Note: A point vector data set named "NORMALS" will be used for normal vectors.

          @param file_name   The name of the output file
          @param description This is included as a file comment at the start of the file.
          @return 0 if everything worked, and non-zero otherwise */
      int write_xml_vtk(std::string file_name, std::string description) { 
        /* Check that we have data */
        if (num_points() <= 0) {
          std::cout << "ERROR(write_xml_vtk): No points!" << std::endl;
          return 1;
        }
        if (num_cells() <= 0) {
          std::cout << "ERROR(write_xml_vtk): No cells!" << std::endl;
          return 2;
        }
        /* Looks like we have data.  Let's open our file */
        std::ofstream outStream;
        outStream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (outStream.is_open()) {
          outStream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_xml_vtk): Could not open file!" << std::endl;
          return 3;
        }
        outStream << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
        outStream << "<!-- " << description << " -->" << std::endl;
        outStream << "  <UnstructuredGrid>" << std::endl;
        outStream << "    <Piece NumberOfPoints='" << num_points() << "' NumberOfCells='" << num_cells() << "'>" << std::endl;
        if ((pdata_sdat.size() > 0) || (pdata_vdat.size() > 0)) {
          bool firstName;
          outStream << "      <PointData Scalars='";
          firstName = true;
          if ( !(pdata_sdat.empty())) {
            for (auto& kv : pdata_sdat) {
              outStream << (firstName ? "" : " ") << kv.first;
              firstName = false;
            }
            outStream << "'";
          }
          firstName = true;
          outStream << " ";
          if ( !(pdata_vdat.empty())) {
            if (pdata_vdat.contains("NORMALS"))
              outStream << "Normals='" << "NORMALS" << "' ";
            outStream << "Vectors='";
            for (auto& kv : pdata_vdat) {
              if (kv.first != "NORMALS") {
                outStream << (firstName ? "" : " ") << kv.first;
                firstName = false;
              }
            }
            outStream << "'";
          }
          outStream << ">" << std::endl;
          if (pdata_sdat.size() > 0) {
            for (auto& kv : pdata_sdat) {
              kv.second.resize(num_points()); // Just in case we don't have data for all points. ;)
              outStream << "        <DataArray Name='" << kv.first << "' type='Float64' format='ascii' NumberOfComponents='1'>" << std::endl;
              outStream << "          ";
              for (const auto& dv : kv.second) 
                outStream << std::setprecision(10) << dv << " ";
              outStream << std::endl;
              outStream << "        </DataArray>" << std::endl;
            }
          }
          if (pdata_vdat.size() > 0) {
            for (auto& kv : pdata_vdat) {
              kv.second.resize(num_points()); // Just in case we don't have data for all points. ;)
              outStream << "        <DataArray Name='" << kv.first << "' type='Float64' format='ascii' NumberOfComponents='3'>" << std::endl;
              for (const auto& dv : kv.second) 
                outStream << "          " << std::setprecision(10) << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
              outStream << "        </DataArray>" << std::endl;
            }
          }
          outStream << "      </PointData>" << std::endl;
        }
        outStream << "      <Points>" << std::endl;
        outStream << "        <DataArray Name='Points' type='Float64' format='ascii' NumberOfComponents='3'>" << std::endl;
        for (const auto& pnt : pnt_idx_to_pnt) 
          outStream << "          " << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
        outStream << "        </DataArray>" << std::endl;
        outStream << "      </Points>" << std::endl;
        outStream << "      <Cells>" << std::endl;
        outStream << "        <DataArray type='Int32' Name='connectivity' format='ascii'>" << std::endl;
        for(auto& poly: cell_lst) {
          outStream << "          " ;
          for(auto& vert: poly)
            outStream << vert << " ";
          outStream << std::endl;
        }
        outStream << "        </DataArray>" << std::endl;
        outStream << "        <DataArray type='Int32' Name='offsets'      format='ascii'>" << std::endl;
        outStream << "          ";
        std::vector<int>::size_type j = 0;
        for(auto& poly: cell_lst) {
          j += poly.size();
          outStream << j << " ";
        }
        outStream << std::endl;
        outStream << "        </DataArray>" << std::endl;
        outStream << "        <DataArray type='Int8' Name='types'      format='ascii'>" << std::endl;
        outStream << "          ";
        for(auto& poly: cell_lst)
          outStream << cell_type_to_vtk_type(req_pt_cnt_to_cell_type(poly.size())) << " ";
        outStream << std::endl;
        outStream << "        </DataArray>" << std::endl;
        outStream << "      </Cells>" << std::endl;
        outStream << "    </Piece>" << std::endl;
        outStream << "  </UnstructuredGrid>" << std::endl;
        outStream << "</VTKFile>" << std::endl;
        
        /* Final newline */
        outStream << std::endl;
        outStream.close();
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to a legacy VTK UNSTRUCTUREDGRID file.

          Note: A point vector data set named "NORMALS" will be used for normal vectors.
          Note: A point vector data set named "COLORS" will be used as COLOR_SCALARS.  Components should be in [0.0, 1.0].

          @param file_name   The name of the output file
          @param description This is the file description.
          @return 0 if everything worked, and non-zero otherwise */
        int write_legacy_vtk(std::string file_name, std::string description) { 
          /* Check that we have data */
          if (num_points() <= 0) {
            std::cout << "ERROR(write_legacy_vtk): No points!" << std::endl;
            return 1;
          }
          if (num_cells() <= 0) {
            std::cout << "ERROR(write_legacy_vtk): No cells!" << std::endl;
            return 2;
          }
          /* Looks like we have data.  Let's open our file */
          std::ofstream outStream;
          outStream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
          if (outStream.is_open()) {
            outStream.imbue(std::locale::classic());
          } else {
            std::cout << "ERROR(write_legacy_vtk): Could not open file!" << std::endl;
            return 3;
          }
          /* Dump the header */
          outStream << "# vtk DataFile Version 3.0" << std::endl;
          outStream << description << std::endl;
          outStream << "ASCII" << std::endl;
          outStream << "DATASET UNSTRUCTURED_GRID" << std::endl;
          /* Dump the points */
          outStream << "POINTS " << num_points() << " double" << std::endl;
          for (const auto& pnt : pnt_idx_to_pnt) {
            outStream << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
          }
          /* Dump the cell data */
          std::vector<int>::size_type total_cells_ints = 0;
          for(auto& cell: cell_lst) 
            total_cells_ints += (1+cell.size());
          outStream << "CELLS " << num_cells() << " " << total_cells_ints << std::endl;
          for(auto& poly: cell_lst) {
            outStream << poly.size() << " ";
            for(auto& vert: poly) {
              outStream << vert << " ";
            }
            outStream << std::endl;
          }
          outStream << "CELL_TYPES " << num_cells() << std::endl;
          for(auto& poly: cell_lst)
            outStream << cell_type_to_vtk_type(req_pt_cnt_to_cell_type(poly.size())) << std::endl;
          /* Dump point scalar data */
          if ((pdata_sdat.size() > 0) || (pdata_vdat.size() > 0)) {
            outStream << "POINT_DATA " << num_points() << std::endl;
            if (pdata_sdat.size() > 0) {
              for (auto& kv : pdata_sdat) {
                kv.second.resize(num_points()); // Just in case we don't have data for all points. ;)
                outStream << "SCALARS " << kv.first << " double 1" << std::endl;
                outStream << "LOOKUP_TABLE default" << std::endl;
                for (const auto& dv : kv.second) {
                  outStream << std::setprecision(10) << dv << std::endl;
                }
              }
            }
            /* Dump point vector data */
            if (pdata_vdat.size() > 0) {
              for (auto& kv : pdata_vdat) {
                kv.second.resize(num_points()); // Just in case we don't have data for all points. ;)
                if ("NORMALS" == kv.first) {
                  outStream << "NORMALS " << kv.first << " double" << std::endl;
                } else if ("COLORS" == kv.first) {
                  outStream << "COLOR_SCALARS " << kv.first << " 3" << std::endl;
                } else {
                  outStream << "VECTORS " << kv.first << " double" << std::endl;
                }
                for (const auto& dv : kv.second) 
                  outStream << std::setprecision(10) << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
              }
            }
          }
          /* Final newline */
          outStream << std::endl;
          outStream.close();
          return 0;
        }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump object to STDOUT. 

          @param max_num_print Maximum number of points/cells to print.  Use 0 to print all points/cells. */
      void dump_cplx(int max_num_print) const {
        std::cout << "Meta Data" << std::endl;
        std::cout << "  Points ............. " << num_points() << std::endl;
        std::cout << "  Scalar Data Sets ... " << pdata_sdat.size() << std::endl;
        std::cout << "  Vector Data Sets ... " << pdata_vdat.size() << std::endl;
        std::cout << "  Cells .............. " << num_cells() << std::endl;
        print_all_points(max_num_print);
        print_all_cells(max_num_print);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to a PLY file

          Note: If one of the point vector data sets is named "NORMALS", then it will be used for normal vectors.
          Note: If one of the point vector data sets is named "COLORS", then it will be used for color data.  Components should be in [0.0, 1.0].

          @param file_name   The name of the output file
          @param description For legacy files, this is the file description.  For XML files this is included as a file comment at the start of the file.
          @return 0 if everything worked, and non-zero otherwise */
      int write_ply(std::string file_name, std::string description) { 
        /* Check that we have data */
        if (num_points() <= 0) {
          std::cout << "ERROR(write_ply): No points!" << std::endl;
          return 1;
        }
        if (num_cells() <= 0) {
          std::cout << "ERROR(write_ply): No cells!" << std::endl;
          return 2;
        }
        for(auto& poly: cell_lst) {
          cell_type_t cell_type = req_pt_cnt_to_cell_type(poly.size());
          if ( !((cell_type == cell_type_t::TRIANGLE) || (cell_type == cell_type_t::QUAD))) {
            std::cout << "ERROR(write_ply): Cells must all be 2D (triangles or quads)!" << std::endl;
            return 2;
          }
        }
        /* Looks like we have data.  Let's open our file */
        std::ofstream outStream;
        outStream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (outStream.is_open()) {
          outStream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_ply): Could not open file!" << std::endl;
          return 3;
        }
        bool have_color_data = pdata_vdat.contains("COLORS");
        bool have_normal_data = pdata_vdat.contains("NORMALS");
        outStream << "ply" << std::endl;
        outStream << "format ascii 1.0" << std::endl;
        outStream << "comment software: Mitch Richling's MR_rect_tree package" << std::endl;
        outStream << "comment note: " << description << std::endl;
        outStream << "element vertex " << num_points() << std::endl;
        outStream << "property float x" << std::endl;
        outStream << "property float y" << std::endl;
        outStream << "property float z" << std::endl;
        if (have_color_data) {
          outStream << "property uchar red" << std::endl;
          outStream << "property uchar green" << std::endl;
          outStream << "property uchar blue" << std::endl;
        }
        if (have_normal_data) {
          outStream << "property float nx" << std::endl;
          outStream << "property float ny" << std::endl;
          outStream << "property float nz" << std::endl;
        }
        outStream << "element face " << num_cells() << std::endl; // May need to be adjusted if cells are not triangles..
        outStream << "property list uchar int vertex_index" << std::endl;
        outStream << "end_header" << std::endl;
        // Dump Vertex Data
        for (int i=0; i<num_points(); i++) {
          pnt_t pnt = pnt_idx_to_pnt[i];
          outStream << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2];
          if (have_color_data) {
            vdat_t clr = pdata_vdat["COLORS"][i];
            outStream << " " << static_cast<int>(255*clr[0]) << " " << static_cast<int>(255*clr[1]) << " " << static_cast<int>(255*clr[2]);
          }
          if (have_normal_data) {
            vdat_t nml = pdata_vdat["NORMALS"][i];
            double nmll = std::sqrt(nml[0]*nml[0]+nml[1]*nml[1]+nml[2]*nml[2]);
            if (nmll < eps)
              nmll = 1.0;
            outStream << " " << std::setprecision(10) << (nml[0]/nmll) << " " << (nml[1]/nmll) << " " << (nml[2]/nmll);
          }
          outStream << std::endl;
        }
        // Dump Cells
        for(auto& poly: cell_lst) {
          outStream << poly.size() << " ";
          for(auto& vert: poly) {
            outStream << vert << " ";
          }
          outStream << std::endl;
        }
        /* Final newline */
        outStream << std::endl;
        outStream.close();
        return 0;
      }
      //@}

  };

  typedef MR_cell_cplx<false, false, false, false, false, 1.0e-3> MRccF3;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  1.0e-3> MRccT3;

  typedef MR_cell_cplx<false, false, false, false, false, 1.0e-5> MRccF5;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  1.0e-5> MRccT5;

  typedef MR_cell_cplx<false, false, false, false, false, 1.0e-9> MRccF9;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  1.0e-9> MRccT9;

}

#define MJR_INCLUDE_MR_cell_cplx
#endif
