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

set(TEST_SOURCES_IO
        iges_reader_writer_test.cc
        io_converter_test.cc
        irit_reader_writer_test.cc
        vtk_writer_test.cc
        xml_reader_test.cc
        xml_writer_test.cc)

set(IGES_READ ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/read.iges)
set(IGES_READ_2 ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/read2.iges)
set(PATH_TO_IRIT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/test.itd)
set(PATH_TO_XML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/io/test_files/spline_tank.xml)
configure_file(config_irit.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_irit.h)
configure_file(config_xml.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_xml.h)
configure_file(config_iges.in.h ${CMAKE_CURRENT_BINARY_DIR}/config_iges.h)
