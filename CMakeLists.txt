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


set(TARGETS_EXAMPLE "MR_cell_cplx_test_mpl" "MR_cell_cplx_test_wedge" 
                    "MR_cell_cplx_test_surf_wave_hole" 
                    "hello_world" 
                    "MR_rect_tree_test_surf_wave" 
                    "MR_rect_tree_test_surf_param" 
                    "MR_rect_tree_test_surf_param_deg"
                    "MR_rect_tree_test_surf_impl" 
                    "MR_rect_tree_test_curve2d_impl" 
                    "MR_rect_tree_test_surf_edge" 
                    "MR_rect_tree_test_surf_corner"
                    "MR_rect_tree_test_flat1"
                    "MR_rect_tree_test_cplx"
                  )
set(TARGETS_TEST "check_cell_hexahedron" "check_cell_pyramid" "check_cell_quad" "check_cell_triangle" "check_cell_segment"
                 "geomi_seg_isect_type" 
               )

# Construct list of targets we can build
set(TARGETS_ALL ${TARGETS_EXAMPLE} ${TARGETS_TEST})

list(REMOVE_DUPLICATES TARGETS_ALL)

if(TARGETS_TEST)
  if(Boost_FOUND)
    enable_testing()
  else()
    list(REMOVE_ITEM TARGETS_ALL ${TARGETS_TEST})
    message(WARNING "Unit tests will NOT be built!")
  endif()
endif()

# Add a target for each example we can build
foreach(CURTGT IN LISTS TARGETS_ALL)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${CURTGT}.cpp")
    add_executable(${CURTGT} "examples/${CURTGT}.cpp")
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
  DEPENDS ${TARGETS_TEST}
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
