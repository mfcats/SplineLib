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

cmake_minimum_required(VERSION 3.6)

project(SplineLibExamples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

#set(SOURCES ex1.cc)
#set(OPTIMIZER testcase_optimizer.cc)
set(IGA iga.cc)

#find_package(SplineLib REQUIRED)
find_package(CSiga REQUIRED)

#add_executable(ex1 ${SOURCES})
#add_executable(TestcaseOptimizer ${OPTIMIZER})
add_executable(iga ${IGA})

#target_link_libraries(ex1 SplineLib::splinelibspl)
#target_link_libraries(TestcaseOptimizer SplineLib::splinelibspl SplineLib::splinelibutil SplineLib::splinelibio)
target_link_libraries(iga CSiga)

