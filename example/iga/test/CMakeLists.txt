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
        basis_function_handler_test.cc
        bdf_handler_test.cc
        connectivity_handler_test.cc
        element_generator_test.cc
        element_integral_calculator_test.cc
        element_integration_point_test.cc
        element_test.cc
        integration_point_test.cc
        integration_rule_test.cc
        linear_equation_assembler_test.cc
        mapping_handler_test.cc
        solution_vtk_writer_examples.cc
        solution_vtk_writer_test.cc)

find_package(GMock REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(CSigaTests ${TEST_SOURCES})
target_link_libraries(CSigaTests CSiga GMock::Main)

target_compile_definitions(CSigaTests PRIVATE GTEST_HAS_PTHREAD=0)
