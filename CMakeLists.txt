######################################################################################################################################################
######################################################################################################################################################
cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

set(CMAKE_CXX_EXTENSIONS FALSE)

######################################################################################################################################################
# Proejct meta data

project(MRTree CXX)

set(PROJECT_VERSION_DATE  2024-06-24)
set(PROJECT_VERSION_ID    v0.0.30)
mark_as_advanced(PROJECT_VERSION_DATE PROJECT_VERSION_ID)

######################################################################################################################################################
# Here are the options we use to modify the build -- these options are passed to cmake like so: cmake -DO_BTEST=YES or -DO_BTEST=NO

OPTION(O_DOXYGEN "Include support for doxygen"         ON )
OPTION(O_BTEST   "Include support for boost.test"      ON )

######################################################################################################################################################
# Find tools & libraries used to build primary/secondary/tertiary targets

if(O_DOXYGEN)
  find_package(Doxygen)
endif()

if(O_BTEST)
find_package(Boost COMPONENTS unit_test_framework)
endif()

######################################################################################################################################################


set(TARGETS_EXAMPLE "hello_world" 
                    "surface_with_normals" "surface_plot_edge" "surface_branch_glue" "surface_plot_corner"
                    "parametric_surface_with_defects" "performance_with_large_surface" "trefoil" 
                    "parametric_curve_3d"
                    "implicit_surface" 
                    "implicit_curve_2d" 
                    "vector_field_3d"
                    "curve_plot"
                    "complex_magnitude_surface"
                    "ear_surface"
                  )
set(TARGETS_UTEST "check_cell_hexahedron" "check_cell_pyramid" "check_cell_quad" "check_cell_triangle" "check_cell_segment"
                 "geomi_pnt_line_distance"
                 "geomi_seg_isect_type"
                 "geomr_pnt_pln_distance" "geomr_pnt_tri_distance" "geomr_pnt_line_distance"
                 "tree_basics_7b1" "tree_basics_7b2" "tree_basics_7b3" "tree_basics_7b4" "tree_basics_7b5"
                 "tree_basics_15b1" "tree_basics_15b3"
                 "tree_corners"
                 "zzz_tbd_tree_children" "zzz_tbd_tree_neighbors"
               )
set(TARGETS_FTEST "MR_cell_cplx_test_mpl" "MR_cell_cplx_test_wedge" "MR_cell_cplx_test_surf_wave_hole" 
                  #"MR_rect_tree_test_surf_corner" 
                  "flat_test_tree_01"
                  "nan_solver"
               )

# Construct list of targets we can build
set(TARGETS_ALL ${TARGETS_EXAMPLE} ${TARGETS_UTEST} ${TARGETS_FTEST})

list(REMOVE_DUPLICATES TARGETS_ALL)

if(TARGETS_UTEST)
  if(Boost_FOUND)
    enable_testing()
  else()
    list(REMOVE_ITEM TARGETS_ALL ${TARGETS_UTEST})
    message(WARNING "Unit tests will NOT be built!")
  endif()
endif()

# Add a target for each example we can build
foreach(CURTGT IN LISTS TARGETS_ALL)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${CURTGT}.cpp")
    add_executable(${CURTGT} "examples/${CURTGT}.cpp")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ftests/${CURTGT}.cpp")
    add_executable(${CURTGT} "ftests/${CURTGT}.cpp")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    add_executable(${CURTGT} "utests/${CURTGT}.cpp")
  else()
    message("Warning: Unable to find source for target ${CURTGT}!")
    continue()
  endif()
  
  target_include_directories(${CURTGT} PRIVATE lib)
  target_compile_features(${CURTGT} PUBLIC cxx_std_23)
  target_include_directories(${CURTGT} PRIVATE ${CMAKE_BINARY_DIR})

  # Compiler specific stuff
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${CURTGT} PUBLIC -Wall -Wconversion -Wno-unknown-pragmas -Wextra -Wno-deprecated-copy)
    target_compile_options(${CURTGT} PUBLIC -O4)
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    target_compile_options(${CURTGT} PUBLIC -Wall -Wconversion -Wno-unknown-pragmas -Wextra -Wno-sign-conversion)
    target_compile_options(${CURTGT} PUBLIC -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_compile_options(${CURTGT} PUBLIC -Wall)
    target_compile_options(${CURTGT} PUBLIC -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    #message("Warning: MSVC support is currently experimental")
  endif()

  # Unit Test Specific Stuff
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    target_link_libraries(${CURTGT} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(NAME "${CURTGT}" COMMAND "${CURTGT}")
  endif()

endforeach(CURTGT)

######################################################################################################################################################

add_custom_target(test-build 
  DEPENDS ${TARGETS_UTEST}
  COMMENT "Building Unit Tests"
)

add_custom_target(examples
  DEPENDS ${TARGETS_EXAMPLE}
  COMMENT "Building Examples"
)

add_custom_target(clean-data
  COMMAND  rm -f *.vtu *.vtp *.vtk *.ply
  COMMENT "Cleaning VTK & PLY files"
)

######################################################################################################################################################
# Generate an include file with various methdata about the build.

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.hpp.cmake ${CMAKE_BINARY_DIR}/config.hpp)

######################################################################################################################################################

# Add targets for documentation
if(Doxygen_FOUND)
  foreach(DOXINPUT IN ITEMS "lib" "examples")
    message("Info: Generateing doxygen target for ${DOXINPUT}")
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}/Doxyfile)
    add_custom_target(doc-${DOXINPUT}
      COMMAND ${DOXYGEN_EXECUTABLE} > dox.out
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}
      COMMENT "Generating documentation with Doxygen"
      VERBATIM)
  endforeach(DOXINPUT)
else()
  message("Warning: Doxygen not found.  No documentation targets!")
endif()
