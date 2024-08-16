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
OPTION(O_MRASTER "Include support for MRaster"         ON )

######################################################################################################################################################
# Find tools & libraries used to build targets

if(O_DOXYGEN)
  find_package(Doxygen)
endif()

if(O_BTEST)
  find_package(Boost COMPONENTS unit_test_framework)
endif()

if(O_MRASTER)
  find_file (MRASTER_PATH "ramCanvasTpl.hpp" PATHS "mraster"  "../mraster" PATH_SUFFIXES "lib")
  if(NOT MRASTER_PATH STREQUAL "MRASTER_PATH-NOTFOUND")
    cmake_path(GET MRASTER_PATH PARENT_PATH MRASTER_PATH)
    set(MRASTER_FOUND "YES")
  endif()
endif()

######################################################################################################################################################

set(TARGETS_EXAMPLE "hello_world" 
                    "surface_with_normals" 
                    "surface_plot_edge" "surface_branch_glue" "surface_plot_corner" "surface_plot_step" "surface_plot_annular_edge"
                    "parametric_surface_with_defects" "performance_with_large_surface" "trefoil" 
                    "parametric_curve_3d"
                    "implicit_surface" 
                    "implicit_curve_2d" 
                    "vector_field_3d"
                    "curve_plot"
                    "complex_magnitude_surface"
                    "ear_surface"
                    "holy_wave_surf"
                  )
set(TARGETS_UTEST "check_cell_hexahedron" "check_cell_pyramid" "check_cell_quad" "check_cell_triangle" "check_cell_segment"
                  "geomi_pnt_line_distance"
                  "geomi_seg_isect_type"
                  "geomr_pnt_pln_distance" "geomr_pnt_tri_distance" "geomr_pnt_line_distance"
                  "tree_basics_7b1" "tree_basics_7b2" "tree_basics_7b3" "tree_basics_7b4" "tree_basics_7b5"
                  "tree_basics_15b1" "tree_basics_15b3"
                  "tree_corners"
                  "tree_children" "tree_neighbors" "zzz_tbd_tree_exist_neighbors"
                )
set(TARGETS_FTEST "flat_test_tree_01"
                  "nan_solver"
                  "segment_folder" "triangle_folder"
                  "rect_fix_dup" "rect_fix_nan"
                )
set(TARGETS_MRASTER  )

# Construct list of targets we can build
set(TARGETS_ALL ${TARGETS_EXAMPLE} ${TARGETS_UTEST} ${TARGETS_FTEST})

list(REMOVE_DUPLICATES TARGETS_ALL)

# If we found boost test, then enable it. If not, then remove targets requiring it.
if(TARGETS_UTEST)
  if(NOT MRASTER_FOUND)
    list(REMOVE_ITEM TARGETS_ALL ${TARGETS_MRASTER})
    message(WARNING "Examples requiring MRaster NOT be built!")
  endif()
endif()

# If we found didn't find MRaster, then remove targets requiring it.
if(TARGETS_MRASTER)
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

  # targets that require MRaster
  if(MRASTER_FOUND)
    list (FIND TARGETS_MRASTER ${CURTGT} TMPR)
    if(${TMPR} GREATER -1)
      target_include_directories(${CURTGT} PRIVATE ${MRASTER_PATH})
    endif()
  endif()

endforeach(CURTGT)

######################################################################################################################################################

add_custom_target(ftest-build 
  DEPENDS ${TARGETS_FTEST}
  COMMENT "Building Functional Tests"
)

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
