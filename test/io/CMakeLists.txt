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

set(TEST_SOURCES
        ${TEST_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/iges_reader_writer_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/io_converter_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/irit_reader_writer_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk_writer_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/xml_reader_test.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/xml_writer_test.cc
        PARENT_SCOPE)

set(IGES_READ ${CMAKE_CURRENT_SOURCE_DIR}/test_files/read.iges PARENT_SCOPE)
set(IGES_READ_2 ${CMAKE_CURRENT_SOURCE_DIR}/test_files/read2.iges PARENT_SCOPE)
set(PATH_TO_IRIT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test_files/test.itd PARENT_SCOPE)
set(PATH_TO_XML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test_files/spline_tank.xml PARENT_SCOPE)
configure_file(config_irit.in.h ${PROJECT_BINARY_DIR}/test/config_irit.h)
configure_file(config_xml.in.h ${PROJECT_BINARY_DIR}/test/config_xml.h)
configure_file(config_iges.in.h ${PROJECT_BINARY_DIR}/test/config_iges.h)
