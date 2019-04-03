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

add_library(splinelibutil INTERFACE)
target_include_directories(splinelibutil INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

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
        element.h
        element_generator.h
        named_type.h
        numeric_operations.h
        numeric_settings.h
        multi_index_handler.h
        random.h
        string_operations.h
        system_operations.h
        vector_utils.h
        DESTINATION "${include_install_dir}")
