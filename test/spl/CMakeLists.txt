# Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University
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

set(
        TEST_SOURCES
        ${TEST_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/b_spline_2d_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/b_spline_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/control_point_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/degree_elevation_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/degree_reduction_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/knot_insertion_and_removal_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/knot_insertion_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/knot_removal_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/nurbs_2d_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/nurbs_3d_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/nurbs_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/parameter_space_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/physical_space_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/projection_curve_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/projection_surface_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/random_b_spline_generator_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/random_nurbs_generator_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/spline_subdivision_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/square_generator_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/surface_generator_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/weighted_physical_space_test.cc
        PARENT_SCOPE
)

set(
        MAX_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/surface_generator_test
        PARENT_SCOPE
)
