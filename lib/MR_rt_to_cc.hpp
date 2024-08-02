// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_rt_to_cc.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Implimentation for the MR_rt_to_cc class.@EOL
 @keywords  tree cell complex
 @std       C++23
 @see       MR_rect_tree.hpp, MR_cell_cplx.hpp
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
#ifndef MJR_INCLUDE_MR_rt_to_cc

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <tuple>                                                         /* STL tuples              C++11    */
#include <vector>                                                        /* STL vector              C++11    */
#include <string>                                                        /* C++ strings             C++11    */
#include <variant>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put everything in the mjr namespace
namespace mjr {
  template <class rt_t, class cc_t>
  /** @brief Tessellates a MR_rect_tree object, and places the result into an MR_cell_cplx object.

      One might think of this class as a fancy pseudo-constructor or proto-factory for MR_cell_cplx objects.  Structurally it's simply a templated collection
      of types and static methods.

      @tparam rt_t The type of the attached MR_rect_tree object
      @tparam cc_t The type of the attached MR_cell_cplx object */
  class MR_rt_to_cc {

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Describe point source. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Specify a source space for a data index. */
      enum class tree_val_src_t { DOMAIN,    //!< The domain space.
                                  RANGE,     //!< The range space.
                                  CONSTANT   //!< A pseudo-source that returns a constant.
                                };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold an integer or double */
      typedef std::variant<int, typename cc_t::uft_t> iorf_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold a description of how to extract a scalar value from a tree object */
      typedef std::tuple<tree_val_src_t, iorf_t> tree_scl_val_desc_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A list of tree_scl_val_desc_t objects.  */
      typedef std::vector<tree_scl_val_desc_t> tree_scl_val_desc_lst_t;
      //@}

    private:
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Utility Functions. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A list of tree_scl_val_desc_t objects.  */
      static void create_dataset_to_point_mapping(const rt_t& rtree, cc_t& ccplx, const tree_scl_val_desc_lst_t& rt_dil) {

        typename cc_t::data_idx_lst_t cc_data_idx_lst(3);
        for(int i=0; i<3; ++i)
          if(get<0>(rt_dil[i]) == tree_val_src_t::DOMAIN)
            cc_data_idx_lst[i] = get<int>(get<1>(rt_dil[i]));
          else if(get<0>(rt_dil[i]) == tree_val_src_t::RANGE)
            cc_data_idx_lst[i] = get<int>(get<1>(rt_dil[i])) + rtree.domain_dimension;
          else if(get<0>(rt_dil[i]) == tree_val_src_t::CONSTANT)
            cc_data_idx_lst[i] = get<typename cc_t::uft_t>(get<1>(rt_dil[i]));
        ccplx.create_dataset_to_point_mapping(cc_data_idx_lst);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt coordinates, extract point/scalar/vector data, and add point/data to cc
          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param diti                 The point coordinate in rtree */
      static typename cc_t::pnt_idx_t add_point_and_data_from_tree(cc_t&        ccplx,
                                                                   const rt_t&  rtree,
                                                                   rt_t::diti_t diti) {
        return add_point_and_data_from_data(ccplx, rtree, rtree.diti_to_drpt(diti), rtree.get_sample(diti));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt coordinates, extract point/scalar/vector data, and add point/data to cc
          @param ccplx    The MR_cell_cplx to populate with geometry
          @param rtree    The MR_rect_tree with source data
          @param dom_pnt  Domain point
          @param rng_pnt  Range point */
      static typename cc_t::pnt_idx_t add_point_and_data_from_data(cc_t&                 ccplx,
                                                                   const rt_t&           rtree,
                                                                   typename rt_t::drpt_t dom_pnt,
                                                                   typename rt_t::rrpt_t rng_pnt) {
        typename cc_t::pnt_data_t pd;
        if constexpr (rtree.domain_dimension == 1)
          pd.push_back(dom_pnt);
        else
          for(auto v: dom_pnt)
            pd.push_back(v);
        if constexpr (rtree.range_dimension == 1)
          pd.push_back(rng_pnt);
        else
          for(auto v: rng_pnt)
            pd.push_back(v);
        return ccplx.add_point(pd);
      }
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Poly data construction */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate attached MR_cell_cplx object from data in attached MR_rect_tree object.

          construct_geometry_fans(), unlike the other geometry construction methods, is capable of "healing" some broken edges -- edges with one good point
          and one NaN point.  It uses the origional sampling function to solve along the edge to produce the longest non-NaN edge possible.  Then it
          constructs cells using the new piont(s).  If func is nullptr, then no edge healing is preformed.  This feature works for segments & triangles only.

          @verbatim
             | Geom       | Dom Dim | Out Dim | Result             |
             |------------+---------+---------+--------------------|
             | POINTS     |     1-3 |       1 | All Cell Points    |
             |------------+---------+---------+--------------------|
             | RECTANGLES |     1-3 |       0 | Cell Corner Points |
             | RECTANGLES |     2-3 |       1 | Cell Edges         |
             | RECTANGLES |       2 |       2 | 2D Rectangles      |
             | RECTANGLES |       3 |       2 | Cell Faces         |
             | RECTANGLES |       3 |       3 | Solid Hexahedra    |
             |------------+---------+---------+--------------------|
             | FANS       |       2 |       1 | Triangle Edges     |
             | FANS       |       2 |       2 | Triangles          |
             | FANS       |       3 |       1 | Pyramid Edges      | TODO: Not currently implemented
             | FANS       |       3 |       2 | Pyramid Faces      | TODO: Not currently implemented
             | FANS       |       3 |       3 | Solid Pyramids     | TODO: Not currently implemented
          @endverbatim

          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param cells                List of cells to output from rtree
          @param output_dimension     Parts of cells to output
          @param point_src            Point sources
          @param func                 The function was used to sample the tree */
      static int construct_geometry_fans(cc_t&                         ccplx,
                                         const rt_t&                   rtree,
                                         typename rt_t::diti_list_t    cells,
                                         int                           output_dimension,
                                         tree_scl_val_desc_lst_t       point_src,
                                         typename rt_t::rsfunc_t       func = nullptr
                                   ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        if (rtree.domain_dimension == 1) {
          for(auto& cell: cells) {
            typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(ccplx, rtree, cell);
            typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell);
            typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[0]);
            typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[1]);
            if (func) { // We have a func, so we can "heal" broken edges.
              if (ctr_pnti < 0) { // Center: Broken. Left: 
                if(cn0_pnti >= 0) { // Center: Broken.  Left: Good.
                  typename cc_t::pnt_idx_t np = nan_edge_solver(ccplx, rtree, cn0_pnti, corners[0], cell, func);
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cn0_pnti, np}, output_dimension);
                }
                if(cn1_pnti >= 0) { // Center: Broken.  Right: Good.
                  typename cc_t::pnt_idx_t np = nan_edge_solver(ccplx, rtree, cn1_pnti, corners[1], cell, func);
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {np, cn1_pnti}, output_dimension);
                }
              } else {             // Center: Good.
                if(cn0_pnti < 0) { // Center: Good.  Left: Broken.
                  typename cc_t::pnt_idx_t np = nan_edge_solver(ccplx, rtree, ctr_pnti, cell, corners[0], func);
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {np, ctr_pnti}, output_dimension);
                } else {           // Center: Good.  Left: Good.
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
                }
                if(cn1_pnti < 0) { // Center: Good.  Right: Broken.
                  typename cc_t::pnt_idx_t np = nan_edge_solver(ccplx, rtree, ctr_pnti, cell, corners[1], func);
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {ctr_pnti, np}, output_dimension);
                } else {           // Center: Good.  Left: Good.
                  ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
                }
              }
            } else {
              ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
              ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
            }
          }
        } else if (rtree.domain_dimension == 2) {
          for(auto& cell: cells) {
            if (func) { // We have a func, so we can "heal" broken edges.
              for(int i=0; i<2; i++) {
                for(int j=-1; j<2; j+=2) {
                  std::vector<typename rt_t::diti_list_t> triangles;
                  typename rt_t::diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                  if (nbrs.size() > 1) {
                    for(auto n: nbrs) {
                      typename rt_t::diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                      if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                        triangles.push_back({corners[1], corners[0], cell});
                      else
                        triangles.push_back({corners[0], corners[1], cell});
                    }
                  } else {
                    typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                    if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                      triangles.push_back({corners[1], corners[0], cell});
                    else
                      triangles.push_back({corners[0], corners[1], cell});
                  }
                  for(auto triangle: triangles) {
                    std::array<typename cc_t::pnt_idx_t, 3> tpnts {add_point_and_data_from_tree(ccplx, rtree, triangle[0]),
                                                                   add_point_and_data_from_tree(ccplx, rtree, triangle[1]),
                                                                   add_point_and_data_from_tree(ccplx, rtree, triangle[2])};
                    int num_bad = static_cast<int>(std::count_if(tpnts.begin(), tpnts.end(), [](typename cc_t::pnt_idx_t i) { return i<0; }));
                    if (num_bad == 0) {
                      ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {tpnts[0], tpnts[1], tpnts[2]}, output_dimension);
                    } else if ((num_bad == 1) || (num_bad == 2)) {
                      // Rotate points so we only have two cases to think about...
                      std::array<int, 3> p {0, 1, 2};
                      if ( ((tpnts[1] < 0) && (num_bad == 1)) || ((tpnts[1] >= 0) && (num_bad == 2)) )
                        p = {1, 2, 0};
                      else if ( ((tpnts[2] < 0) && (num_bad == 1)) || ((tpnts[2] >= 0) && (num_bad == 2)) )
                        p = {2, 0, 1};
                      // Solve for edge 0-1 & 0-2
                      if (num_bad == 1) {
                        typename cc_t::pnt_idx_t np1 = nan_edge_solver(ccplx, rtree, tpnts[p[1]], triangle[p[1]], triangle[p[0]], func);
                        typename cc_t::pnt_idx_t np2 = nan_edge_solver(ccplx, rtree, tpnts[p[2]], triangle[p[2]], triangle[p[0]], func);
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {np1, tpnts[p[1]], tpnts[p[2]]}, output_dimension);
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {tpnts[p[2]], np2, np1}, output_dimension);
                      } else {
                        typename cc_t::pnt_idx_t np1 = nan_edge_solver(ccplx, rtree, tpnts[p[0]], triangle[p[0]], triangle[p[1]], func);
                        typename cc_t::pnt_idx_t np2 = nan_edge_solver(ccplx, rtree, tpnts[p[0]], triangle[p[0]], triangle[p[2]], func);
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {tpnts[p[0]], np1, np2}, output_dimension);
                      }
                    }
                  }
                }
              }
            } else { // We don't a func, so we can can't "heal" broken edges.  This is much faster. ;)
              typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(ccplx, rtree, cell);
              if (ctr_pnti >= 0) { // Center point was good, let's try and make some triangles...
                for(int i=0; i<2; i++) {
                  for(int j=-1; j<2; j+=2) {
                    typename rt_t::diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                    if (nbrs.size() > 1) {
                      for(auto n: nbrs) {
                        typename rt_t::diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                        typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[0]);
                        typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[1]);
                        if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                          std::swap(cn0_pnti, cn1_pnti);
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                      }
                    } else {
                      typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                      typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[0]);
                      typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(ccplx, rtree, corners[1]);
                      if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                        std::swap(cn0_pnti, cn1_pnti);
                      ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                    }
                  }
                }
              }
            }

          }
        } else if (rtree.domain_dimension == 3) {
          std::cout << "ERROR: construct_geometry_fans: domain_dimension==3 not supported!" << std::endl;
          return 1;
          //  MJR TODO NOTE <2024-08-02T09:43:13-0500> construct_geometry_fans: Implement
        } else { //if (rtree.domain_dimension > 3) {
          std::cout << "ERROR: construct_geometry_fans: output_dimension>3 not supported for output_dimension>0!" << std::endl;
          return 1;
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry_fans(cc_t&                         ccplx,
                                         const rt_t&                   rtree,
                                         int                           output_dimension,
                                         tree_scl_val_desc_lst_t       point_src,
                                         typename rt_t::rsfunc_t                func = nullptr
                                   ) {
        return construct_geometry_fans(ccplx,
                                       rtree,
                                       rtree.get_leaf_cells(rtree.ccc_get_top_cell()),
                                       output_dimension,
                                       point_src,
                                       func);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate attached MR_cell_cplx object from data in attached MR_rect_tree object.

          Only 0D vertex cells are produced.  While vertex cells may be obtained by setting the output_dimension to zero and calling construct_geometry_fan or
          construct_geometry_rects, this method is much faster.  This method also provides the option of only outputting centers -- no corners.

          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param cells                List of cells to output from rtree
          @param point_src            Point sources
          @param output_centers       Create vertexes for cell  centers
          @param output_corners       Create vertexes for cell corners*/
      static int construct_geometry_points(cc_t&                         ccplx,
                                           const rt_t&                   rtree,
                                           typename rt_t::diti_list_t    cells,
                                           tree_scl_val_desc_lst_t       point_src,
                                           bool                          output_centers,
                                           bool                          output_corners
                                          ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        if (output_centers && output_corners) {
          for(auto& cell: cells)
            for(auto& vert: rtree.ccc_get_vertexes(cell))
              ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(ccplx, rtree, vert)});
        } else if (output_centers) {
          for(auto& cell: cells)
            ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(ccplx, rtree, cell)});
        } else if (output_corners) {
          for(auto& cell: cells)
            for(auto& vert: rtree.ccc_get_corners(cell))
              ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(ccplx, rtree, vert)});
        } else {
          std::cout << "WARNING: construct_geometry_points: Both output_centers & output_corners are FALSE.  No geometry created!" << std::endl;
          return 1;
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      static int construct_geometry_rects(cc_t&                         ccplx,
                                          const rt_t&                   rtree,
                                          typename rt_t::diti_list_t    cells,
                                          int                           output_dimension,
                                          tree_scl_val_desc_lst_t       point_src
                                         ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        for(auto& cell: cells) {
          std::vector<typename cc_t::pnt_idx_t> cnr_pti;
          typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell);
          for(auto& corner: corners) {
            typename cc_t::pnt_idx_t pnti = add_point_and_data_from_tree(ccplx, rtree, corner);
            cnr_pti.push_back(pnti);
          }
          if (rtree.domain_dimension == 1) {
            ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cnr_pti[0], cnr_pti[1]}, output_dimension);
          } else if (rtree.domain_dimension == 2) {
            ccplx.add_cell(cc_t::cell_type_t::QUAD, {cnr_pti[0], cnr_pti[1], cnr_pti[3], cnr_pti[2]}, output_dimension);
          } else { // if(rtree.domain_dimension == 3) {
            ccplx.add_cell(cc_t::cell_type_t::HEXAHEDRON,
                           {cnr_pti[0], cnr_pti[1], cnr_pti[3], cnr_pti[2],
                            cnr_pti[4], cnr_pti[5], cnr_pti[7], cnr_pti[6]},
                           output_dimension);
          }
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry_rects(cc_t&                    ccplx,
                                          const rt_t&              rtree,
                                          int                      output_dimension,
                                          tree_scl_val_desc_lst_t  point_src
                                         ) {
        return construct_geometry_rects(ccplx,
                                  rtree,
                                  rtree.get_leaf_cells(rtree.ccc_get_top_cell()),
                                  output_dimension,
                                  point_src);
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Mathematical Tools */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given an edge with one good point and one NaN point, find the longest segment from the good point toward the NaN point.

          Note we normally use this function when we detect a NaN in a geometric point (i.e. the things with a pnt_idx_t).  This solver solves until
          the return of func has no NaNs.  Those two criteria might not be the same thing, but it's OK.

          @param ccplx                  The MR_cell_cplx to populate with geometry
          @param rtree                  The MR_rect_tree with source data
          @param good_point_ccplx_index Good point index in the ccplx object
          @param good_point_rtree_index Good point index in the rtree object
          @param sick_point_rtree_index Bad point index in the rtree object
          @param func                   The function to use for the solver */
      static typename cc_t::pnt_idx_t nan_edge_solver(cc_t&                    ccplx,
                                                      const rt_t&              rtree,
                                                      typename cc_t::pnt_idx_t good_point_ccplx_index,
                                                      typename rt_t::diti_t    good_point_rtree_index,
                                                      typename rt_t::diti_t    sick_point_rtree_index,
                                                      typename rt_t::rsfunc_t  func,
                                                      typename cc_t::uft_t     solver_epsilon=cc_t::epsilon/100
                                                     ) {
        // Solver cache.  Clear it if we have a different rtree object from last time.
        static std::unordered_map<typename rt_t::diti_t, std::unordered_map<typename rt_t::diti_t, typename cc_t::pnt_idx_t>> nan_solver_cache;
        static const rt_t* rtree_cache = nullptr;
        if (rtree_cache != &rtree) {
          nan_solver_cache.clear();
          rtree_cache = &rtree;
        }
        // Check to see if we solved this one before
        if (nan_solver_cache.contains(sick_point_rtree_index))
          if (nan_solver_cache[sick_point_rtree_index].contains(good_point_rtree_index))
            return  nan_solver_cache[sick_point_rtree_index][good_point_rtree_index];
        // Apparently we need to solve this one as it's not in the case
        typename rt_t::drpt_t good_point_drpt = rtree.diti_to_drpt(good_point_rtree_index);
        typename rt_t::drpt_t sick_point_drpt = rtree.diti_to_drpt(sick_point_rtree_index);
        typename rt_t::rrpt_t good_point_rrpt = rtree.get_sample(good_point_rtree_index);
        typename rt_t::drpt_t init_point_drpt = good_point_drpt;
        while ( (rtree.drpt_distance_inf(good_point_drpt, sick_point_drpt) > solver_epsilon) ) {
          typename rt_t::drpt_t md_point_drpt = rtree.drpt_midpoint(good_point_drpt, sick_point_drpt);
          typename rt_t::rrpt_t y = func(md_point_drpt);
          if (rtree.rrpt_is_nan(y)) {
            sick_point_drpt = md_point_drpt;
          } else {
            good_point_drpt = md_point_drpt;
            good_point_rrpt = y;
          }
        }
        // Figure out what to return, add it to the cache, and return.
        typename cc_t::pnt_idx_t ret;
        if (rtree.drpt_distance_inf(good_point_drpt, init_point_drpt) < (ccplx.epsilon)) // Use ccplx here!!!
          ret = good_point_ccplx_index;
        else
          ret = add_point_and_data_from_data(ccplx, rtree, good_point_drpt, good_point_rrpt);
        nan_solver_cache[sick_point_rtree_index][good_point_rtree_index] = ret;
        return ret;
      }
      //@}

  };
}

#define MJR_INCLUDE_MR_rt_to_cc
#endif
