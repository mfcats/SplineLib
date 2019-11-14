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

cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(SplineLib VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(namespace "${PROJECT_NAME}::")
set(targets_export_name "${PROJECT_NAME}Targets")

# Installation
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
include(CMakePackageConfigHelpers)

set(include_install_dir "include")
set(library_install_dir "lib")
set(executable_install_dir "lib")
set(config_install_dir "lib/cmake/${PROJECT_NAME}")

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Compile using clang compiler!")
    set(OPTIMIZATION_FLAGS -O3 -march=native)
    set(RUNTIME_CHECKS_DEBUG "")
    set(PARALLELIZATION_FLAGS "")
    set(WARNING_FLAGS -Wall -Wextra -Wpedantic -Wunused -Wno-missing-braces) # -Werror)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    message(STATUS "Compile using gcc compiler!")
    set(OPTIMIZATION_FLAGS -O3 -march=native)
    set(RUNTIME_CHECKS_DEBUG "")
    set(PARALLELIZATION_FLAGS "")
    set(WARNING_FLAGS -Wall -Wextra -Wpedantic -Wunused) # -Werror)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(STATUS "Compile using Intel compiler!")
    set(OPTIMIZATION_FLAGS -fast)
    set(RUNTIME_CHECKS_DEBUG -fp-trap=all -fp-trap-all=all -ftrapuv)
    set(PARALLELIZATION_FLAGS -vec -parallel)
    set(WARNING_FLAGS -Wall -Wcheck -Weffc++) # -Werror)
else()
    message(ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
set(CMAKE_DEBUG_POSTFIX "d")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

# TODO(all): Is this how it should to be as DEBUG seems not to be set automatically with CMAKE_BUILD_TYPE=Debug?
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif ()

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(test)

# Package configuration
configure_package_config_file("cmake/Config.cmake.in" "${project_config}" INSTALL_DESTINATION "${config_install_dir}")
write_basic_package_version_file("${version_config}" COMPATIBILITY SameMajorVersion)

install(
        FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}"
)
install(
        EXPORT "${targets_export_name}"
        NAMESPACE "${namespace}"
        DESTINATION "${config_install_dir}"
)
