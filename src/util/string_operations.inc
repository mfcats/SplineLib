/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

#include <string>
#include <vector>

template<typename TYPE>
TYPE ConvertStringToNumber(std::string string) {
  int sign = 1;
  if (string[0] == '-') {
    sign = -1;
    string.erase(0, 1);
  }
  int const found_dot = string.find_first_of('.');
  int exponent = string.find_first_of("Ee");
  auto const npos = static_cast<int>(std::string::npos);
  int const string_length = string.length();
  int end_of_number;
  if (exponent == npos) {
    end_of_number = (string_length - 1);
  } else {
    end_of_number = (exponent - 1);
  }
  double result = 0, factor;
  if (found_dot != npos) {
    factor = pow(10, found_dot - 1);
  } else if (exponent == npos) {
    factor = pow(10, string_length - 1);
  } else {
    factor = pow(10, exponent - 1);
  }
  for (int i = 0; i <= end_of_number; ++i) {
    if (i != found_dot) {
      result += (std::stoi(string.substr(i, 1)) * factor);
      factor /= 10;
    }
  }
  if (exponent != npos) {
    if (string[exponent + 1] == '-') {
      factor = -1;
      ++exponent;
    } else if (string[exponent + 1] == '+') {
      ++exponent;
    }
    int const potency = std::stoi(string.substr(exponent + 1, string_length - exponent));
    if (factor != -1) {
      result *= pow(10, potency);
    } else {
      result /= pow(10, potency);
    }
  }
  result *= sign;
  return static_cast<TYPE>(result);
}

template<typename TYPE>
std::vector<TYPE> ConvertStringVectorToNumberVector(std::vector<std::string> const &string_vector) {
  std::vector<TYPE> converted;
  converted.reserve(string_vector.size());
  for (std::string const &string : string_vector) converted.emplace_back(ConvertStringToNumber<double>(string));
  return converted;
}

template<typename TYPE>
std::vector<TYPE> ConvertDelimitedStringToNumberVector(std::string string) {
  std::vector<TYPE> vector;
  while (!string.empty()) {
    int const found_comma = string.find_first_of(',');
    int const found_semicolon = string.find_first_of(';');
    if (((found_comma < found_semicolon) || (found_semicolon == -1)) && (found_comma > 0)) {
      vector.push_back(ConvertStringToNumber<TYPE>(TrimSpacesAndSquareBrackets(string.substr(0, found_comma))));
      string.erase(0, found_comma + 1);
    } else if (((found_semicolon < found_comma) || (found_comma == -1)) && (found_semicolon > 0)) {
      vector.push_back(ConvertStringToNumber<TYPE>(TrimSpacesAndSquareBrackets(string.substr(0, found_semicolon))));
      string.erase(0, found_semicolon + 1);
    } else {
      string.erase(0, 1);
    }
  }
  return vector;
}
