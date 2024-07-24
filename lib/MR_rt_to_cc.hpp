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
      /** @name Describe geometric structure. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Geometric Structure */
      enum class cell_structure_t { RECTANGLES,   //!< hyper-rectangles made up of a cell's corners.
                                    FANS,         //!< hyper-triangles containing the center point and shared corner points of neighboring cells
                                  };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Describe tree value sources. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Specify a source space for a data index. */
      enum class tree_val_src_t { DOMAIN,    //!< The domain space.
                                  RANGE,     //!< The range space.
                                  CONSTANT   //!< A pseudo-source that returns a constant.
                                };

      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold a description of how to extract a scalar value from a tree object */
      typedef std::variant<int, typename cc_t::pnt_crd_t> iorf_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold a description of how to extract a scalar value from a tree object */
      typedef std::tuple<typename cc_t::pdata_name_t, tree_val_src_t, iorf_t> tree_scl_val_desc_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A list of tree_scl_val_desc_t objects.  */
      typedef std::vector<tree_scl_val_desc_t> tree_scl_val_desc_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold a description of how to extract a 3-vector value from a tree object. */
      typedef std::tuple<typename cc_t::pdata_name_t, tree_val_src_t, iorf_t, tree_val_src_t, iorf_t, tree_val_src_t, iorf_t> tree_vec_val_desc_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A list of tree_vec_val_desc_t objects. */
      typedef std::vector<tree_vec_val_desc_t> tree_vec_val_desc_lst_t;
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Describe tree value sources. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a tree_val_src to a string */
      static std::string tree_val_src_to_string(tree_val_src_t tree_val_src) { 
        switch(tree_val_src) {
          case tree_val_src_t::DOMAIN:   return("DOMAIN"); break;
          case tree_val_src_t::RANGE:    return("RANGE");  break;
          case tree_val_src_t::CONSTANT: return("CONSTANT");   break;
        }
        return("ERROR"); // Never get here if the switch above is correct..
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a tree_scl_val_desc_t to a string */
      static std::string tree_scl_val_desc_to_string(tree_scl_val_desc_t scl_desc) { 
        if (std::get<2>(scl_desc).index() == 0)
          return (std::get<0>(scl_desc) + ":" + tree_val_src_to_string(get<1>(scl_desc)) + "/" + std::to_string(get<int>(get<2>(scl_desc))));
        else
          return (std::get<0>(scl_desc) + ":" + tree_val_src_to_string(get<1>(scl_desc)) + "/" + std::to_string(get<typename cc_t::pnt_crd_t>(get<2>(scl_desc))));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Validate a tree_scl_val_desc_t object
          @param rtree    The MR_rect_tree with source data
          @param scl_desc Object to valudate */
      static bool invalid_tree_scl_val_desc(const rt_t& rtree, tree_scl_val_desc_t scl_desc) { 
        if (std::get<1>(scl_desc) == tree_val_src_t::CONSTANT) {
          if (std::get<2>(scl_desc).index() != 1) {
            std::cout << "ERROR: Invalid val_desc (const must be double): " << tree_scl_val_desc_to_string(scl_desc) << std::endl;
            return true;
          }
        } else {
          if (std::get<2>(scl_desc).index() != 0) {
            std::cout << "ERROR: Invalid val_desc (idx must be int): " << tree_scl_val_desc_to_string(scl_desc) << std::endl;
            return true;
          }
          if (std::get<int>(std::get<2>(scl_desc)) < 0) {
            std::cout << "ERROR: Invalid val_desc (idx too small): " << tree_scl_val_desc_to_string(scl_desc) << std::endl;
            return true;
          }
          if (std::get<1>(scl_desc) == tree_val_src_t::DOMAIN) {
            if (std::get<int>(std::get<2>(scl_desc)) >= rtree.domain_dimension) {
              std::cout << "ERROR: Invalid val_desc (idx too large): " << tree_scl_val_desc_to_string(scl_desc) << std::endl;
              return true;
            }
          } else {
            if (std::get<int>(std::get<2>(scl_desc)) >= rtree.range_dimension) {
              std::cout << "ERROR: Invalid val_desc (idx too large): " << tree_scl_val_desc_to_string(scl_desc) << std::endl;
              return true;
            }
          }
        }
        return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Validate a tree_vec_val_desc_t object
          @param rtree    The MR_rect_tree with source data
          @param vec_desc Object to valudate */
      static bool invalid_tree_vec_val_desc(const rt_t& rtree, tree_vec_val_desc_t vec_desc) { 
        return (invalid_tree_scl_val_desc(rtree, {std::get<0>(vec_desc), std::get<1>(vec_desc), std::get<2>(vec_desc)}) ||
                invalid_tree_scl_val_desc(rtree, {std::get<0>(vec_desc), std::get<3>(vec_desc), std::get<4>(vec_desc)}) ||
                invalid_tree_scl_val_desc(rtree, {std::get<0>(vec_desc), std::get<5>(vec_desc), std::get<6>(vec_desc)}));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract data described by a tree_scl_val_desc_t
          @param rtree    The MR_rect_tree with source data
          @param index         The input tree_scl_val_desc_t object
          @param domainValue   Data for the domain value
          @param rangeValue    Data for the range value */
      static typename cc_t::pnt_crd_t get_scalar(const rt_t& rtree, tree_scl_val_desc_t index, typename rt_t::drpt_t& domainValue, typename rt_t::rrpt_t& rangeValue) {
        if (std::get<1>(index) == tree_val_src_t::DOMAIN) {
          return (rtree.dom_at(domainValue, std::get<int>(std::get<2>(index))));
        } else if (std::get<1>(index) == tree_val_src_t::RANGE) {
          return (rtree.rng_at(rangeValue, std::get<int>(std::get<2>(index))));
        } else if (std::get<1>(index) == tree_val_src_t::CONSTANT) {
          return (std::get<typename cc_t::pnt_crd_t>(std::get<2>(index)));
        } else {
          std::cout << "ERROR: Invalid index source" << std::endl;
          return (0.0);
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract data described by a tree_vec_val_desc_t
          @param rtree         The MR_rect_tree with source data
          @param index         The input tree_vec_val_desc_t object
          @param domainValue   Data for the domain value
          @param rangeValue    Data for the range value */
      static typename cc_t::pnt_t get_vector(const rt_t& rtree, tree_vec_val_desc_t index, typename rt_t::drpt_t& domainValue, typename rt_t::rrpt_t& rangeValue) {
        return typename cc_t::pnt_t({get_scalar(rtree, tree_scl_val_desc_t({std::get<0>(index), std::get<1>(index), std::get<2>(index)}), domainValue, rangeValue),
                                     get_scalar(rtree, tree_scl_val_desc_t({std::get<0>(index), std::get<3>(index), std::get<4>(index)}), domainValue, rangeValue),
                                     get_scalar(rtree, tree_scl_val_desc_t({std::get<0>(index), std::get<5>(index), std::get<6>(index)}), domainValue, rangeValue)});
      }


      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt coordinates, extract point/scalar/vector data, and add point/data to cc 
          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param diti                 The point coordinate in rtree
          @param point_src            Point descriptor
          @param scalar_data_src_lst  Scalar data descriptors
          @param vector_data_src_lst  Vector data descriptors */
      static typename cc_t::pnt_idx_t add_point_and_data_from_tree(const rt_t&  rtree, 
                                                                   cc_t&        ccplx,
                                                                   rt_t::diti_t diti,
                                                                   tree_vec_val_desc_t     const& point_src, 
                                                                   tree_scl_val_desc_lst_t const& scalar_data_src_lst, 
                                                                   tree_vec_val_desc_lst_t const& vector_data_src_lst) {
        typename rt_t::drpt_t    dPts = rtree.diti_to_drpt(diti);
        typename rt_t::rrpt_t    rPts = rtree.get_sample(diti);
        typename cc_t::pnt_t     pnt  = get_vector(rtree, point_src, dPts, rPts);
        typename cc_t::pnt_idx_t pnti = ccplx.add_point(pnt);
        if (ccplx.last_point_added_was_new()) { // Logically unnecessary, but speeds things up.
          for (auto s : scalar_data_src_lst) 
            ccplx.add_data_if_new(std::get<0>(s), get_scalar(rtree, s, dPts, rPts));
          for (auto v : vector_data_src_lst)
            ccplx.add_data_if_new(std::get<0>(v), get_vector(rtree, v, dPts, rPts));
        }
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
             | FANS       |       3 |       1 | Pyramid Edges      |
             | FANS       |       2 |       2 | Triangles          |
             | FANS       |       3 |       2 | Pyramid Faces      |
             | FANS       |       3 |       3 | Solid Pyramids     |
          @endverbatim

          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param cells                List of cells to output from rtree
          @param cell_structure       Type of cell data to output
          @param output_dimension     Parts of cells to output
          @param point_src            Point sources
          @param scalar_data_src_lst  List of point data sources (scalars)
          @param vector_data_src_lst  List of point data sources (vectors) */
      static int construct_geometry(cc_t&                         ccplx,
                                    const rt_t&                   rtree, 
                                    typename rt_t::diti_list_t    cells,
                                    cell_structure_t              cell_structure,
                                    int                           output_dimension,
                                    tree_vec_val_desc_t           point_src, 
                                    tree_scl_val_desc_lst_t       scalar_data_src_lst, 
                                    tree_vec_val_desc_lst_t       vector_data_src_lst
                                   ) {
        if (output_dimension < 0) {
          std::cout << "ERROR: construct_geometry: output_dimension < 0!" << std::endl;
          return 1;
        }
        if (invalid_tree_vec_val_desc(rtree, point_src)) {
          return 1;
        }
        if (std::any_of(scalar_data_src_lst.cbegin(), scalar_data_src_lst.cend(), 
                        [rtree](auto s) { return invalid_tree_scl_val_desc(rtree, s); })) {
          return 1;
        }
        if (std::any_of(vector_data_src_lst.cbegin(), vector_data_src_lst.cend(), 
                        [rtree](auto v) { return invalid_tree_vec_val_desc(rtree, v); })) {
          return 1;
        }
        ccplx.clear();
        if (output_dimension == 0) {
          for(auto& cell: cells) {
            typename rt_t::diti_list_t verts = (cell_structure == cell_structure_t::FANS ? rtree.ccc_get_vertexes(cell) : rtree.ccc_get_corners(cell));
            for(auto& vert: verts) {
              typename cc_t::pnt_idx_t pnti = add_point_and_data_from_tree(rtree, ccplx, vert, point_src, scalar_data_src_lst, vector_data_src_lst);
              ccplx.add_cell(cc_t::cell_type_t::POINT, {pnti});
            }
          }
        } else { // if(output_dimension > 0) {
          if (rtree.domain_dimension == 3) {
            std::cout << "ERROR: construct_geometry: output_dimension >3 not supported for output_dimension>0!" << std::endl;
            return 1;
          }
          if (cell_structure == cell_structure_t::FANS) {
            if (rtree.domain_dimension == 1) {
              for(auto& cell: cells) {
                typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(rtree, ccplx, cell, point_src, scalar_data_src_lst, vector_data_src_lst);
                typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell);
                typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0], point_src, scalar_data_src_lst, vector_data_src_lst);
                typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1], point_src, scalar_data_src_lst, vector_data_src_lst);
                ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
                ccplx.add_cell(cc_t::cell_type_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
              }
            } else if (rtree.domain_dimension == 2) {
              for(auto& cell: cells) {
                typename cc_t::pnt_idx_t ctr_pnti = add_point_and_data_from_tree(rtree, ccplx, cell, point_src, scalar_data_src_lst, vector_data_src_lst);
                for(int i=0; i<2; i++) {
                  for(int j=-1; j<2; j+=2) {
                    typename rt_t::diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                    if (nbrs.size() > 1) {
                      for(auto n: nbrs) {
                        typename rt_t::diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                        typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0], point_src, scalar_data_src_lst, vector_data_src_lst);
                        typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1], point_src, scalar_data_src_lst, vector_data_src_lst);
                        if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                          ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn1_pnti, cn0_pnti, ctr_pnti}, output_dimension);
                        else
                          ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                      }
                    } else {
                      typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                      typename cc_t::pnt_idx_t cn0_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[0], point_src, scalar_data_src_lst, vector_data_src_lst);
                      typename cc_t::pnt_idx_t cn1_pnti = add_point_and_data_from_tree(rtree, ccplx, corners[1], point_src, scalar_data_src_lst, vector_data_src_lst);
                      if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn1_pnti, cn0_pnti, ctr_pnti}, output_dimension);
                      else
                        ccplx.add_cell(cc_t::cell_type_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                    }
                  }
                }
              }
            } else { // if (rtree.domain_dimension == 3) {
              if (output_dimension == 3) {
                //  MJR TODO NOTE construct_geometry: ADD CODE
              } else if (output_dimension == 2) {
                //  MJR TODO NOTE construct_geometry: ADD CODE
              } else { // (output_dimension == 1)
                //  MJR TODO NOTE construct_geometry: ADD CODE
              }
            }          
          } else { // if (cell_structure == cell_structure_t::RECTANGLES) {
            for(auto& cell: cells) {
              std::vector<typename cc_t::pnt_idx_t> cnr_pti;
              typename rt_t::diti_list_t corners = rtree.ccc_get_corners(cell);
              for(auto& corner: corners) {
                typename cc_t::pnt_idx_t pnti = add_point_and_data_from_tree(rtree, ccplx, corner, point_src, scalar_data_src_lst, vector_data_src_lst);
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
          }
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry(cc_t&                           ccplx,
                                    const rt_t&                     rtree, 
                                    cell_structure_t                cell_structure,
                                    int                             output_dimension,
                                    tree_vec_val_desc_t             point_src, 
                                    tree_scl_val_desc_lst_t         scalar_data_src_lst, 
                                    tree_vec_val_desc_lst_t         vector_data_src_lst
                                   ) {
        return construct_geometry(ccplx,
                                  rtree,
                                  rtree.get_leaf_cells(rtree.ccc_get_top_cell()),
                                  cell_structure,
                                  output_dimension,
                                  point_src, 
                                  scalar_data_src_lst, 
                                  vector_data_src_lst);
      }
      //@}

  };
}

#define MJR_INCLUDE_MR_rt_to_cc
#endif
 

