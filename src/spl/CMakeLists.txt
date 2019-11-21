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
        SOURCES
        control_point.cc
        surface_generator.cc
        square_generator.cc
)

add_library(splinelibspl SHARED ${SOURCES})
target_include_directories(splinelibspl PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)
target_link_libraries(splinelibspl PUBLIC splinelibbaf splinelibutil)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(splinelibspl PUBLIC ${RUNTIME_CHECKS_DEBUG} ${WARNING_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(splinelibspl PUBLIC ${OPTIMIZATION_FLAGS} ${WARNING_FLAGS} ${PARALLELIZATION_FLAGS})
endif()

install(
        TARGETS splinelibspl
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(
        FILES
        b_spline.h
        b_spline_generator.h
        control_point.h
        control_point.inc
        nurbs.h
        nurbs_generator.h
        parameter_space.h
        physical_space.h
        projection.h
        random_b_spline_generator.h
        random_nurbs_generator.h
        random_spline_utils.h
        spline.h
        spline_generator.h
        square_generator.h
        surface_generator.h
        weighted_physical_space.h
        DESTINATION "${include_install_dir}/src/spl"
)
