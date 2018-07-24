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
        element.cc
        element_generator.cc
        element_integration_point.cc)

add_library(splinelibelm SHARED ${SOURCES})
target_include_directories(splinelibelm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(splinelibelm splinelibbaf)

install(
        TARGETS splinelibelm
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)


install(FILES
        element.h
        element_generator.h
        element_integration_point.h
        DESTINATION "${include_install_dir}")

