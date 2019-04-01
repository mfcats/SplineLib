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

add_subdirectory(executables)

set(SOURCES
        converter_log.cc
        iges_reader.cc
        iges_writer.cc
        irit_reader.cc
        irit_writer.cc
        reader.cc
        vtk_writer.cc
        writer.cc
        xml_reader.cc
        xml_writer.cc)

add_library(splinelibio SHARED ${SOURCES})

target_include_directories(splinelibio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(splinelibio splinelibspl pugixmllib)

install(
        TARGETS splinelibio
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
        converter_log.h
        iges_reader.h
        iges_writer.h
        io_converter.h
        irit_reader.h
        irit_reader_utils.h
        irit_writer.h
        irit_writer_utils.h
        reader.h
        vtk_writer.h
        vtk_writer_utils.h
        xml_reader.h
        xml_reader_utils.h
        xml_writer.h
        xml_writer_utils.h
        writer.h
        DESTINATION "${include_install_dir}")
