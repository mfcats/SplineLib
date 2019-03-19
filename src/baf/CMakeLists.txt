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

set(SOURCES
        b_spline_basis_function.cc
		basis_function.cc
        basis_function_factory.cc
        control_point.cc
        knot_vector.cc
		zero_degree_b_spline_basis_function.cc)

add_library(splinelibbaf SHARED ${SOURCES})
target_include_directories(splinelibbaf PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(splinelibbaf splinelibutil)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_options(splinelibbaf PUBLIC ${RUNTIME_CHECKS_DEBUG} ${WARNING_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
	target_compile_options(splinelibbaf PUBLIC ${OPTIMIZATION_FLAGS} ${WARNING_FLAGS} ${PARALLELIZATION_FLAGS})
endif ()

install(
        TARGETS splinelibbaf
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
        b_spline_basis_function.h
        basis_function.h
        basis_function_factory.h
        control_point.h
        knot_vector.h
        zero_degree_b_spline_basis_function.h
        DESTINATION "${include_install_dir}")
