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
        iges_reader.cc
        iges_writer.cc
        irit_reader.cc
        irit_writer.cc
        vtk_writer.cc
        xml_reader.cc
        xml_writer.cc)

add_library(splinelibio SHARED ${SOURCES})
add_library(ioconversion INTERFACE)
add_executable(iges2irit iges_to_irit.cc converter_log.cc)
add_executable(iges2vtk iges_to_vtk.cc converter_log.cc)
add_executable(iges2xml iges_to_xml.cc converter_log.cc)
add_executable(irit2iges irit_to_iges.cc converter_log.cc)
add_executable(irit2vtk irit_to_vtk.cc converter_log.cc)
add_executable(irit2xml irit_to_xml.cc converter_log.cc)
add_executable(xml2iges xml_to_iges.cc converter_log.cc)
add_executable(xml2irit xml_to_irit.cc converter_log.cc)
add_executable(xml2vtk xml_to_vtk.cc converter_log.cc)

set_target_properties(iges2irit iges2vtk iges2xml irit2iges irit2vtk irit2xml xml2iges xml2irit xml2vtk PROPERTIES ENABLE_EXPORTS 1)

target_include_directories(splinelibio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(splinelibio splinelibspl pugixmllib)
target_link_libraries(iges2irit splinelibio)
target_link_libraries(iges2vtk splinelibio)
target_link_libraries(iges2xml splinelibio)
target_link_libraries(irit2iges splinelibio)
target_link_libraries(irit2vtk splinelibio)
target_link_libraries(irit2xml splinelibio)
target_link_libraries(xml2iges splinelibio)
target_link_libraries(xml2irit splinelibio)
target_link_libraries(xml2vtk splinelibio)
target_link_libraries(ioconversion INTERFACE iges2irit iges2vtk iges2xml irit2iges irit2vtk irit2xml xml2iges xml2irit xml2vtk)

install(
        TARGETS splinelibio ioconversion iges2irit iges2vtk iges2xml irit2iges irit2vtk irit2xml xml2iges xml2irit xml2vtk
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "${library_install_dir}"
        ARCHIVE DESTINATION "${library_install_dir}"
        RUNTIME DESTINATION "${executable_install_dir}"
        INCLUDES DESTINATION "${include_install_dir}"
)

install(FILES
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
