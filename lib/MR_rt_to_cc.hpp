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
      static typename cc_t::pnt_idx_t add_point_and_data_from_tree(const rt_t&  rtree, 
                                                                   cc_t&        ccplx,
                                                                   rt_t::diti_t diti) {



        typename rt_t::drta_t     dPts = rtree.diti_to_drta(diti);
        typename rt_t::rrta_t     rPts = rtree.get_sample_rrta(diti);
        typename cc_t::pnt_data_t pd;
        pd.insert(std::end(pd), std::begin(dPts), std::end(dPts));
        pd.insert(std::end(pd), std::begin(rPts), std::end(rPts));
        typename cc_t::pnt_idx_t  pnti = ccplx.add_point(pd);
        return (pnti);
      }
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Poly data construction */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate attached MR_cell_cplx object from data in attached MR_rect_tree object.

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
          @param point_src            Point sources */
      static int construct_geometry_fans(cc_t&                         ccplx,
                                         const rt_t&                   rtree, 
                                         typename rt_t::diti_list_t    cells,
                                         int                           output_dimension,
                                         tree_scl_val_desc_lst_t       point_src
                                         
                                   ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        if (rtree.domain_dimension == 1) {
          for(auto& cell: cells) {
            typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(rtree, ccplx, cell);
            typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell);
            typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0]);
            typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1]);
            ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
            ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
          }
        } else if (rtree.domain_dimension == 2) {
          for(auto& cell: cells) {
            typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(rtree, ccplx, cell);
            for(int i=0; i<2; i++) {
              for(int j=-1; j<2; j+=2) {
                typename rt_t::diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                if (nbrs.size() > 1) {
                  for(auto n: nbrs) {
                    typename rt_t::diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                    typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0]);
                    typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1]);
                    if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                      ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn1_pnti, cn0_pnti, ctr_pnti}, output_dimension);
                    else
                      ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                  }
                } else {
                  typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                  typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0]);
                  typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1]);
                  if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                    ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn1_pnti, cn0_pnti, ctr_pnti}, output_dimension);
                  else
                    ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                }
              }
            }
          }
        } else if (rtree.domain_dimension == 3) {
          std::cout << "ERROR: construct_geometry_fans: domain_dimension==3 not supported!" << std::endl;
          return 1;
          //  MJR TODO NOTE construct_geometry_fans: Implement
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
                                         tree_scl_val_desc_lst_t        point_src
                                   ) {
        return construct_geometry_fans(ccplx,
                                       rtree,
                                       rtree.get_leaf_cells(rtree.ccc_get_top_cell()),
                                       output_dimension,
                                       point_src);
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
              ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(rtree, ccplx, vert)});
        } else if (output_centers) {
          for(auto& cell: cells)
            ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(rtree, ccplx, cell)});
        } else if (output_corners) {
          for(auto& cell: cells)
            for(auto& vert: rtree.ccc_get_corners(cell)) 
              ccplx.add_cell(cc_t::cell_type_t::POINT, {add_point_and_data_from_tree(rtree, ccplx, vert)});
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
            typename cc_t::pnt_idx_t pnti = add_point_and_data_from_tree(rtree, ccplx, corner);
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

  };
}

#define MJR_INCLUDE_MR_rt_to_cc
#endif
 


 // This is how we can figure out where the domain data is stored in the MR_cell_cplx object.  This could then be exported
 // to the MR_cell_cplx object so that it can do various computations on the mesh...
        // // Figure out how to construct MR_rect_tree domain points from MR_cell_cplx point data sets.
        // std::array<int, rtree.domain_dimension> dom_indexes;
        // dom_indexes.fill(-1);
        // for(int i=0; i<static_cast<int>(scalar_data_src_lst.size()); ++i) 
        //   if (std::get<1>(scalar_data_src_lst[i]) == tree_val_src_t::DOMAIN)
        //     dom_indexes[std::get<int>(std::get<2>(scalar_data_src_lst[i]))] = i;
        // if (std::any_of(dom_indexes.cbegin(), dom_indexes.cend(), [](int i) { return (i<0); }))
        //   std::cout << "WARNING: construct_geometry_fans: All domain variables not captured.  Unable to heal broken edges!" << std::endl;
        // // Traverse the cells and construct geometry
