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

add_library(splinelibio INTERFACE)
target_include_directories(splinelibio INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(splinelibio INTERFACE splinelibspl pugixmllib)

install(
        TARGETS splinelibio
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
