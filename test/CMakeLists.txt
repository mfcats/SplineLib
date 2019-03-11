# Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
#
# This file is part of SplineLib.
#
# SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation version 3 of the License.
#
# SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
# <http://www.gnu.org/licenses/>.
#

set(TEST_SOURCES
        baf/knot_vector_test.cc
        baf/control_point_test.cc
        baf/zero_degree_b_spline_basis_function_test.cc
        baf/b_spline_basis_function_test.cc
        spl/b_spline_test.cc
        spl/b_spline_2d_test.cc
        spl/square_generator_test.cc
        util/multi_index_handler_test.cc
        iga/elm/element_test.cc
        io/iges_reader_writer_test.cc
        iga/itg/integration_rule_test.cc
        iga/itg/integration_point_test.cc
        spl/projection_curve_test.cc
        spl/projection_surface_test.cc
        spl/nurbs_test.cc
        spl/nurbs_2d_test.cc
        spl/nurbs_3d_test.cc
        spl/parameter_space_test.cc
        spl/physical_space_test.cc
        spl/weighted_physical_space_test.cc
        baf/basis_function_factory.cc
        iga/elm/element_integration_point_test.cc
        io/irit_reader_writer_test.cc
        io/vtk_writer_test.cc
        io/xml_reader_test.cc
        io/xml_writer_test.cc
        io/io_converter_test.cc
        iga/connectivity_handler_test.cc
        iga/basis_function_handler_test.cc
        iga/mapping_handler_test.cc
        spl/surface_generator_test.cc
        spl/knot_insertion_test.cc
        spl/nurbs_knot_removal_test.cc
        spl/knot_removal_test.cc
        spl/spline_subdivision_test.cc
        spl/random_b_spline_generator_test.cc
        spl/random_nurbs_generator_test.cc
        iga/element_integral_calculator_test.cc
        iga/linear_equation_assembler_test.cc
        iga/elm/element_generator_test.cc
        iga/solution_vtk_writer_test.cc
        iga/solution_vtk_writer_examples.cc)

set(MAX_SOURCES
  spl/surface_generator_test)

find_package(GMock REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(SplineLibTests ${TEST_SOURCES})
#target_link_libraries(SplineLibTests splinelibbaf splinelibspl splinelibutil splinelibio GMock::Main)
target_link_libraries(SplineLibTests splinelibbaf splinelibiga splinelibspl splinelibutil splinelibio GMock::Main)

add_executable(MaxTest ${MAX_SOURCES})
target_link_libraries(MaxTest splinelibbaf splinelibiga splinelibspl splinelibutil splinelibio GMock::Main)

set(ACCEPTANCE_TEST_SOURCES
        acceptance_tests.cc)
add_executable(AcceptanceTests ${ACCEPTANCE_TEST_SOURCES})
target_link_libraries(AcceptanceTests splinelibbaf splinelibspl splinelibutil)

target_compile_definitions(SplineLibTests PRIVATE GTEST_HAS_PTHREAD=0)
target_compile_definitions(MaxTest PRIVATE GTEST_HAS_PTHREAD=0)

set(IGES_READ ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/read.iges)
set(IGES_READ_2 ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/read2.iges)
set(PATH_TO_IRIT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/test.itd)
set(PATH_TO_XML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/spline_tank.xml)
configure_file(io/config_irit.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_irit.h)
configure_file(io/config_xml.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_xml.h)
configure_file(io/config_iges.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_iges.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
