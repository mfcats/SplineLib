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
        string_operations.cc
        system_operations.cc)

add_library(splinelibutil SHARED ${SOURCES})
target_include_directories(splinelibutil PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

install(
        TARGETS splinelibutil
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
        any_casts.h
        any_casts.inc
        element.h
        element_generator.h
        multi_index_handler.h
        multi_index_handler.inc
        named_type.h
        named_type.inc
        numeric_operations.h
        numeric_operations.inc
        numeric_settings.h
        numeric_settings.inc
        random.h
        random.inc
        stl_container_access.h
        stl_container_access.inc
        string_operations.h
        string_operations.inc
        system_operations.h
        vector_utils.h
        vector_utils.inc
        DESTINATION "${include_install_dir}/src/util")
