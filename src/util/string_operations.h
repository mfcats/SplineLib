/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_STRING_OPERATIONS_H_
#define SRC_UTIL_STRING_OPERATIONS_H_

#include <cmath>
#include <string>
#include <vector>

#include "stl_container_access.h"

namespace splinelib::src::util::string_operations {
bool StartsWith(std::string const &string, std::string const &start_of_string);
bool EndsWith(std::string const &string, std::string const &end_of_string);

std::vector<std::string> SplitStringAtDelimiter(std::string const &string, char delimiter);
std::string TrimSpacesAndSquareBrackets(std::string string);

template<typename TYPE>
TYPE ConvertStringToNumber(std::string string);
template<typename TYPE>
std::vector<TYPE> ConvertStringVectorToNumberVector(std::vector<std::string> const &string_vector);
template<typename TYPE>
std::vector<TYPE> ConvertDelimitedStringToNumberVector(std::string string);

#include "src/util/string_operations.inc"
}  // namespace splinelib::src::util::string_operations

#endif  // SRC_UTIL_STRING_OPERATIONS_H_
