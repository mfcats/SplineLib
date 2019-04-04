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

add_subdirectory(itg)

list(APPEND CMAKE_PREFIX_PATH /usr/local/opt/spack/linux-ubuntu18.04-x86_64/gcc-8/splinelib-github-t7by6dejarqmydvm2xo5zgcqtvkpz3dz)

find_package(Armadillo REQUIRED)
find_package(SplineLib REQUIRED)

add_library(CSiga INTERFACE)
target_include_directories(CSiga INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(CSiga INTERFACE SplineLib::splinelibspl ${ARMADILLO_LIBRARIES})

install(
        TARGETS CSiga
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
        solution_spline.h
        bdf_handler.h
        element_integral_calculator.h
        element_generator.h
        element.h
        linear_equation_assembler.h
        basis_function_handler.h
        connectivity_handler.h
        solution_vtk_writer.h
        mapping_handler.h
        poisson_problem.h
        element_integration_point.h
        DESTINATION "${include_install_dir}")

