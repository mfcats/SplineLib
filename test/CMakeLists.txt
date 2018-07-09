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
        knot_vector_test.cc
        control_point_test.cc
        zero_degree_b_spline_basis_function_test.cc
        b_spline_basis_function_test.cc
        b_spline_test.cc
        2d_b_spline_test.cc
        square_generator_test.cc
        multihandler_test.cc
        element_test.cc
        integration_rule_test.cc
        integration_point_test.cc
        projection_curve_test.cc
        nurbs_test.cc
        nurbs_2_d_test.cc
        nurbs_3_d_test.cc
        parameter_space_test.cc
        physical_space_test.cc
        weighted_physical_space_test.cc)

find_package(GTest REQUIRED)

add_executable(SplineLibTests ${TEST_SOURCES})
target_link_libraries(SplineLibTests splinelibbaf splinelibelm splinelibitg splinelibspl splinelibutil GTest::Main)

target_compile_definitions(SplineLibTests PRIVATE GTEST_HAS_PTHREAD=0)
