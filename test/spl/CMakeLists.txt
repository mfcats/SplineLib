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

set(TEST_SOURCES_SPL
        b_spline_2d_test.cc
        b_spline_test.cc
        knot_insertion_and_removal_test.cc
        knot_insertion_test.cc
        knot_removal_test.cc
        nurbs_2d_test.cc
        nurbs_3d_test.cc
        nurbs_test.cc
        parameter_space_test.cc
        physical_space_test.cc
        projection_curve_test.cc
        projection_surface_test.cc
        random_b_spline_generator_test.cc
        random_nurbs_generator_test.cc
        spline_subdivision_test.cc
        square_generator_test.cc
        surface_generator_test.cc
        weighted_physical_space_test.cc)

set(MAX_SOURCES
        surface_generator_test)
