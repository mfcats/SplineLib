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

add_library(CSigaitg INTERFACE)
target_include_directories(CSigaitg INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

install(
        TARGETS CSigaitg
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
        five_point_gauss_legendre.h
        four_point_gauss_legendre.h
        integration_point.h
        integration_rule.h
        one_point_gauss_legendre.h
        three_point_gauss_legendre.h
        two_point_gauss_legendre.h
        DESTINATION "${include_install_dir}")
