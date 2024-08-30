######################################################################################################################################################
######################################################################################################################################################
cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

######################################################################################################################################################
# Proejct meta data

project(MRPTree 
  VERSION 0.0.0.60
  DESCRIPTION "Generalized bitree/quadtree/octree library"
  HOMEPAGE_URL "https://github.com/richmit/MRPTree/"
  LANGUAGES CXX)

set(MRPTREE_VERSION_DATE  2024-08-20)
mark_as_advanced(MRPTREE_VERSION_DATE)

string(TIMESTAMP MRPTREE_CONFIG_DATE "%Y-%m-%d")
mark_as_advanced(MRPTREE_CONFIG_DATE)

######################################################################################################################################################
# Here are the options we use to modify the build -- these options are passed to cmake like so: cmake -DO_BTEST=YES or -DO_BTEST=NO

OPTION(O_DOXYGEN "Include support for doxygen"         ON )
OPTION(O_BTEST   "Include support for boost.test"      ON )
OPTION(O_MRASTER "Include support for MRaster"         ON )

######################################################################################################################################################
# Find tools & libraries used to build targets

if(O_DOXYGEN)
  set(MRPTREE_OPT_DOXYGEN "1")
  find_package(Doxygen)
  if(Doxygen_FOUND)
    set(MRPTREE_FOUND_DOXYGEN "1")
  else()
    set(MRPTREE_FOUND_DOXYGEN "0")
  endif()
else()
  set(MRPTREE_OPT_DOXYGEN "0")
endif()

if(O_BTEST)
  set(MRPTREE_OPT_BTEST "1")
  find_package(Boost COMPONENTS unit_test_framework)
  if(Boost_FOUND)
    set(MRPTREE_FOUND_BTEST "1")
  else()
    set(MRPTREE_FOUND_BTEST "0")
  endif()
else()
  set(MRPTREE_OPT_BTEST "0")
endif()

if(O_MRASTER)
  set(MRPTREE_OPT_MRASTER "1")
  find_package(MRaster HINTS "./mraster/build/install/share/cmake/MRaster/" "./mraster/build/" "../mraster/build/install/share/cmake/MRaster/" "../mraster/build/")
  message(STATUS "MRaster_FOUND: ${MRaster_FOUND}")
  if(MRaster_FOUND)
    set(MRaster_FOUND "YES")
    set(MRPTREE_FOUND_MRASTER "1")
    set(MRASTER_FOUND_WHAT "PACKAGE")
    message(STATUS "MRaster package found (${MRaster_DIR})")
  else()
    find_file (MRASTER_PATH "CMakeLists.txt" PATHS "mraster"  "../mraster")
    if(NOT MRASTER_PATH STREQUAL "MRASTER_PATH-NOTFOUND")
      cmake_path(GET MRASTER_PATH PARENT_PATH MRASTER_PATH)
      if(EXISTS "${MRASTER_PATH}/build/MRasterLib.cmake")
        set(MRaster_FOUND "YES")
        set(MRPTREE_FOUND_MRASTER "1")
        set(MRASTER_FOUND_WHAT "EXPORT")
        include("${MRASTER_PATH}/build/MRasterLib.cmake")
        message(STATUS "MRaster build package found")
      elseif(EXISTS "${MRASTER_PATH}/lib/ramCanvasTpl.hpp")
        set(MRaster_FOUND "YES")
        set(MRPTREE_FOUND_MRASTER "1")
        set(MRASTER_FOUND_WHAT "SOURCE")
        message(STATUS "MRaster source headers found")
      else()
        set(MRPTREE_FOUND_MRASTER "0")
      endif()
    else()
      set(MRPTREE_OPT_MRASTER "0")
    endif()
  endif()
else()
  set(MRPTREE_OPT_MRASTER "0")
endif()

mark_as_advanced(MRASTER_FOUND_WHAT)
mark_as_advanced(MRPTREE_FOUND_MRASTER MRPTREE_FOUND_BTEST MRPTREE_FOUND_DOXYGEN)
mark_as_advanced(MRPTREE_OPT_MRASTER   MRPTREE_OPT_BOOST   MRPTREE_OPT_DOXYGEN)

######################################################################################################################################################
# Create Targets for examples, functional tests, & unit tests

# Master list of targets
set(TARGETS "hello_world"
            "surface_with_normals"
            "surface_plot_edge" "surface_branch_glue" "surface_plot_corner" "surface_plot_step" "surface_plot_annular_edge"
            "parametric_surface_with_defects" "performance_with_large_surface"
            "parametric_curve_3d"
            "implicit_surface" "implicit_curve_2d"
            "vector_field_3d"
            "curve_plot"
            ### Extra examples
            "ear_surface" "ear_surface_glue"
            "trefoil"
            ### MR_cell_cplx examples (don't use MR_rect_tree or MR_rt_to_cc)
            "holy_wave_surf"
            ### Examples using MRaster
            "hello_MRaster"
            "complex_magnitude_surface"
            ### Unit Tests
            "check_cell_hexahedron" "check_cell_pyramid" "check_cell_quad" "check_cell_triangle" "check_cell_segment"
            "geomi_pnt_line_distance" "geomi_seg_isect_type"
            "geomr_pnt_pln_distance" "geomr_pnt_tri_distance" "geomr_pnt_line_distance"
            "tree_basics_7b1" "tree_basics_7b2" "tree_basics_7b3" "tree_basics_7b4" "tree_basics_7b5"
            "tree_basics_15b1" "tree_basics_15b3"
            "tree_corners"
            "tree_children" "tree_neighbors" "zzz_tbd_tree_exist_neighbors"
            ### Functional Tests
            "flat_test_tree_01"
            "nan_solver"
            "segment_folder" "triangle_folder"
            "rect_fix_dup" "rect_fix_nan"
            ### Config
            "mrptree_config"
           )

# Targets listed here require MRaster to build
set(REQ_MRASTER "hello_MRaster"
                "complex_magnitude_surface"
               )

# Add a target for each one we can build
foreach(CURTGT IN LISTS TARGETS)

  list (FIND REQ_MRASTER ${CURTGT} TMPR)
  if(${TMPR} GREATER -1)
    set(CURTGT_REQ_MRASTER "YES")
  else()
    set(CURTGT_REQ_MRASTER "NO")
  endif()

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${CURTGT}.cpp")
    if(CURTGT_REQ_MRASTER AND NOT MRaster_FOUND)
      message(NOTICE "Target ${CURTGT} will not be built (no MRaster support)")
      continue()
    endif()
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "examples/${CURTGT}.cpp")
    list(APPEND TARGETS_EXAMPLES ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ftests/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "ftests/${CURTGT}.cpp")
    list(APPEND TARGETS_FTEST ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    if(NOT Boost_FOUND)
      message(NOTICE "Target ${CURTGT} will not be built (no boost support)")
      continue()
    endif()
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "utests/${CURTGT}.cpp")
    list(APPEND TARGETS_UTEST ${CURTGT})
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
  if(CURTGT_REQ_MRASTER)
    if(MRASTER_FOUND_WHAT STREQUAL "SOURCE")
      target_include_directories(${CURTGT} PRIVATE "${MRASTER_PATH}/lib")
    else()
      target_link_libraries(${CURTGT} PRIVATE MRaster)
    endif()
  endif()

endforeach(CURTGT)

######################################################################################################################################################

add_custom_target(ftests
  DEPENDS ${TARGETS_FTEST}
  COMMENT "Building Functional Tests"
)

add_custom_target(utests
  DEPENDS ${TARGETS_UTEST}
  COMMENT "Building Unit Tests"
)

add_custom_target(examples
  DEPENDS ${TARGETS_EXAMPLES}
  COMMENT "Building Examples"
)

add_custom_target(clean-data
  COMMAND  rm -f *.vtu *.vtp *.vtk *.ply
  COMMENT "Cleaning VTK & PLY files"
)

if(EXISTS "c:/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-lib/")
  add_custom_target(install-doc-lib
    DEPENDS doc-lib
    COMMAND  rsync -rlt --log-format=%f --stats --delete --delete-excluded --modify-window=2 doc-lib/ /c/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-lib
    COMMENT "Put doc-lib data in web site directory"
  )
endif()

if(EXISTS "c:/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-examples/")
  add_custom_target(install-doc-examples
    DEPENDS doc-examples/
    COMMAND  rsync -rlt --log-format=%f --stats --delete --delete-excluded --modify-window=2 doc-examples/ /c/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-examples
    COMMENT "Put doc-examples data in web site directory"
  )
endif()

######################################################################################################################################################
# Generate an include file with various methdata about the build.

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/mrptree_config.hpp.cmake ${CMAKE_BINARY_DIR}/mrptree_config.hpp)

######################################################################################################################################################

# Add targets for documentation
if(Doxygen_FOUND)
  foreach(DOXINPUT IN ITEMS "lib" "examples")
    message(STATUS "Generateing doxygen target for ${DOXINPUT}")
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}/Doxyfile)
    add_custom_target(doc-${DOXINPUT}
      COMMAND ${DOXYGEN_EXECUTABLE} > dox.out
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}
      COMMENT "Generating ${DOXINPUT} documentation with Doxygen"
      VERBATIM)
  endforeach(DOXINPUT)
else()
  message("Warning: Doxygen not found.  No documentation targets!")
endif()
