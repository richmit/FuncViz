
cc_t::uft_t uc(cc_t::node_data_t xvec) {
  return (xvec[0]-.1);
  //return (1-xvec[0]*xvec[0]-xvec[1]*xvec[1]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t u(tt_t::drpt_t xvec) {
  return (xvec[0] - 0.1);
}

  //ccplx.cull_cells([&ccplx](cc_t::cell_t c){ return ccplx.cell_cross_sdf_boundry(c, uc); });

ccplx.cull_cells([&ccplx](cc_t::cell_t c){ return ccplx.cell_cross_sdf_boundry(c, [](cc_t::node_data_t pd) { return (tc_t::tsdf_to_csdf(u, pd)); }); });






      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Fold triangles that cross over an SDF function.
          @param func           Data function
          @param sdf_index      Index of SDF value in func return
          @param solve_epsilon  Used to detect SDF value near zero */
      void triangle_folder(pdfunc_t func, int sdf_index, uft_t solve_epsilon=epsilon/10) {
        //  MJR TODO NOTE <2024-08-06T12:38:33-0500> triangle_folder: Implement!
        clear_cache_edge_solver_sdf();
        const std::array<int, 3> p00 {0, 1, 2};
        const std::array<int, 3> p10 {1, 2, 0};
        const std::array<int, 3> p20 {2, 0, 1};
        int num_start_cells = num_cells();
        for(int cell_idx=0; cell_idx<num_start_cells; ++cell_idx) {
          if (cell_lst[cell_idx].size() == req_pt_cnt_to_cell_kind(cell_kind_t::TRIANGLE)) {
            auto& cur_cell = cell_lst[cell_idx];
            // Find and count zeros, positive, and negative vertexes
            std::vector<int> zeros;
            std::vector<int> pluss;
            int zero_cnt=0, plus_cnt=0;
            int zero_loc=-1, plus_loc=-1, negv_loc=-1;
            for(int i=0; i<3; i++) {
              node_data_t pd = node_idx_to_node_data[cur_cell[i]];
              if (std::abs(pd[sdf_index]) <= solve_epsilon) {
                zeros[i] = true;
                zero_cnt++;
                zero_loc = i;
              } else {
                zeros[i] = false;
                if (pd[sdf_index] < static_cast<uft_t>(0.0)) {
                  pluss[i] = true;
                  plus_cnt++;
                  plus_loc = i;
                } else {
                  pluss[i] = false;
                  negv_loc = i;
                }
              }
            }
            std::array<int, 3> p;
            if (zero_cnt == 0) { // three triangles or NOP
              if ( (plus_cnt == 1) || (plus_cnt == 2) ) { // three triangles
                // Figure out permutation
                if (plus_cnt == 1) {
                  if (plus_loc == 1)
                    p = p10;
                  else if (plus_loc == 2)  // Rotate 2=>0
                    p = p20;
                  else 
                    p = p00;
                } else {
                  if (negv_loc == 1)
                    p = p10;
                  else if (negv_loc == 2)  // Rotate 2=>0
                    p = p20;
                  else 
                    p = p00;
                }
                // Construct new triangles
                auto orgv0 = cur_cell[0];
                auto orgv1 = cur_cell[1];
                auto orgv2 = cur_cell[2];
                auto newv1 = edge_solver_sdf(orgv0, orgv1, func, sdf_index, solve_epsilon);
                auto newv2 = edge_solver_sdf(orgv0, orgv2, func, sdf_index, solve_epsilon);
                if ((newv1 != orgv0) && (newv1 != orgv1) && (newv2 != orgv0) && (newv2 != orgv2)) { // We got TWO new points
                  cur_cell[1] = newv1; // Modify current triangle in place
                  cur_cell[2] = newv2; // Modify current triangle in place
                  add_cell(cell_kind_t::TRIANGLE, {newv1, orgv1, newv2});  // Add new triangle
                  add_cell(cell_kind_t::TRIANGLE, {orgv1, orgv2, newv2});  // Add new triangle
                }
              }
            } else if (zero_cnt == 1) { // two triangles or NOP
              // Figure out permutation
              if (plus_cnt == 1) { // two triangles
                if (zero_loc == 1)
                  p = p10;
                else if (zero_loc == 2)  // Rotate 2=>0
                  p = p20;
                else 
                  p = p00;
              }
              // Construct new triangles
              auto orgv0 = cur_cell[0];
              auto orgv1 = cur_cell[1];
              auto orgv2 = cur_cell[2];
              auto newv0 = edge_solver_sdf(orgv1, orgv2, func, sdf_index, solve_epsilon);
              if ((newv0 != orgv1) && (newv0 != orgv2)) { // We got new point
                cur_cell[2] = newv0; // Modify current triangle in place
                add_cell(cell_kind_t::TRIANGLE, {orgv0, newv0, orgv2});  // Add new triangle
              }
            }
          }
        }
      }

      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Find the "middle" point of an edge, and add the point to the master point list.

          The geometric midpoint data set is computed from the input points, and this is passed through func (if it is provided).

          @param pnt_idx1      First edge vertex
          @param pnt_idx2      Second edge vertex
          @param func          Mapping function */
      node_idx_t edge_solver_middle(node_idx_t pnt_idx1, node_idx_t pnt_idx2, pdfunc_t func=nullptr) {
        node_data_t pnt_dat1 = node_idx_to_node_data[pnt_idx1];
        node_data_t pnt_dat2 = node_idx_to_node_data[pnt_idx2];
        node_data_t mid_node_data;
        for(decltype(pnt_dat1.size()) i=0; i<pnt_dat1.size(); i++)
          mid_node_data.push_back((pnt_dat1[i] + pnt_dat2[i])/static_cast<uft_t>(2.0));
        if (func)
          mid_node_data = func(mid_node_data);
        return add_node(mid_node_data);
      }



      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add new cells by cutting existing cells with a level plane (a plane orthogonal to one of the axes).

          Return is the number of cells cut.

          UNDER DEVELOPMENT!  Only works on triangle cells!

          New points have interpolated data values.  

          Given a function of point_data_t -> point_data_t.  We could use bisection to solve for a point_data_t such that the constructed point
          from the result is on the level.  We could build a function adapter in MR_rt_to_cc that could provide such a function based on the
          original sample function used by the tree object MR_rt_to_cc created the MR_cell_cplx with.
          

          @param axis_index    Which axis to compare to the level
          @param level         Level to test aginst
          @param close_epsilon Epsilon used to check for "closeness". */
      int cut_on_axis_plane(int axis_index, uft_t level, uft_t close_epsilon=epsilon) {
        int num_cut = 0;
        int num_start_cells = num_cells();
        for(int cell_idx=0; cell_idx<num_start_cells; ++cell_idx) {
          if (cell_lst[cell_idx].size() == 3) {
            std::vector<int> vertex_vs_level_status;
            int stat_plus=0, stat_minus=0, stat_zero=0;
            int loc_plus=0, loc_minus=0, loc_zero=0;
            for(auto pidx: cell_lst[cell_idx]) {
              int tmp = pnt_vs_level(node_idx_to_pnt[pidx], axis_index, level, close_epsilon);
              if (tmp == 0) {
                stat_zero++;
                loc_zero = i;
              } else if (tmp < 0) {
                stat_minus++;
                loc_minus = i;
              } else {
                stat_plus++;
                loc_plus = i;
              }
              vertex_vs_level_status.push_back(tmp);
            }
            if (cell_lst[cell_idx].size() == 3) {
              if ((stat_minus > 0) && (stat_plus > 0)) { // CHOP CHOP
                num_cut++;
                // Figure out which permutation case we need
                int perm_case = 0;
                if (stat_zero == 1)          // rotate to have zero at spot 0
                  perm_case = loc_zero;
                else if (stat_minus == 1)    // rotate to have minus in spot 0
                  perm_case = loc_minus;
                else if (stat_plus == 1)     // rotate to have plus in spot 0
                  perm_case = loc_plus;
                // Compute permutation vector
                std::array<int, 3> p {0, 1, 2};
                if (perm_case == 1)
                  p = {1, 2, 0};
                else if (perm_case == 2)
                  p = {2, 0, 1};
                if (stat_zero == 1) { // Cut into two -- through vertex 0 to opposite side                  
                  uft_t interpolation_constant = intersection distance ratio;
                  new_node_data  = point_data_combination(interpolation_constant, node_idx_to_pnt[cell_list[p[1]]], (1-interpolation_constant), node_idx_to_pnt[cell_list[p[2]]]);
                  new_pnt_idx = add_node(new_node_data);
                  add_cell(cell_list[p[0]], cell_list[p[1]], new_pnt_idx);
                  add_cell(cell_list[p[0]], new_pnt_idx, cell_list[p[2]]);
                  // Zap old cell, or replace existing cell with ^
                } else {              // Cut into three triangles.
                  new_pnt_idx1 = 
                    new_pnt_idx2 =
                    add_cell(cell_list[p[0]], new_pnt_idx1, new_pnt_idx2);
                  add_cell(new_pnt_idx1, cell_list[p[1]], cell_list[p[2]]);
                  add_cell(cell_list[p[2]], new_pnt_idx2, new_pnt_idx1);
                  // Zap old cell, or replace existing cell with ^
                }
              }
            }
          }
        }         
        return num_cut;
      }
