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

add_subdirectory(baf)
add_subdirectory(iga)
add_subdirectory(io)
add_subdirectory(spl)
add_subdirectory(util)

find_package(GMock REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(SplineLibTests ${TEST_SOURCES})
target_link_libraries(SplineLibTests splinelibbaf splinelibiga splinelibspl splinelibutil splinelibio GMock::Main)

add_executable(MaxTest ${MAX_SOURCES})
target_link_libraries(MaxTest splinelibbaf splinelibiga splinelibspl splinelibutil splinelibio GMock::Main)

set(ACCEPTANCE_TEST_SOURCES
        acceptance_tests.cc)
add_executable(AcceptanceTests ${ACCEPTANCE_TEST_SOURCES})
target_link_libraries(AcceptanceTests splinelibbaf splinelibspl splinelibutil)

target_compile_definitions(SplineLibTests PRIVATE GTEST_HAS_PTHREAD=0)
target_compile_definitions(MaxTest PRIVATE GTEST_HAS_PTHREAD=0)

configure_file(io/config_iges.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_iges.h)
configure_file(io/config_irit.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_irit.h)
configure_file(io/config_xml.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_xml.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
