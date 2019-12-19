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
  if (GetFront(string) == '-') {
    sign = -1;
    string.erase(0, 1);
  }
  int const found_dot_index = string.find_first_of('.');
  int exponent_starting_index = string.find_first_of("Ee");
  auto const npos = static_cast<int>(std::string::npos);
  int const string_length = string.length();
  int last_mantissa_index;
  if (exponent_starting_index == npos) {
    last_mantissa_index = (string_length - 1);
  } else {
    last_mantissa_index = (exponent_starting_index - 1);
  }
  double result = 0, factor;
  if (found_dot_index != npos) {
    factor = pow(10, found_dot_index - 1);
  } else if (exponent_starting_index != npos) {
    factor = pow(10, exponent_starting_index - 1);
  } else {
    factor = pow(10, string_length - 1);
  }
  for (int current_digit_index = 0; current_digit_index <= last_mantissa_index; ++current_digit_index) {
    if (current_digit_index != found_dot_index) {
      result += (std::stoi(string.substr(current_digit_index, 1)) * factor);
      factor /= 10;
    }
  }
  if (exponent_starting_index != npos) {
    if (GetValue(string, exponent_starting_index + 1) == '-') {
      factor = -1;
      ++exponent_starting_index;
    } else if (GetValue(string, exponent_starting_index + 1) == '+') {
      ++exponent_starting_index;
    }
    int const potency = std::stoi(string.substr(exponent_starting_index + 1, string_length - exponent_starting_index));
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

// TODO(all): is it always , and ;?
template<typename TYPE>
std::vector<TYPE> ConvertDelimitedStringToNumberVector(std::string string) {
  std::vector<TYPE> vector;
  while (!string.empty()) {
    int const comma_index = string.find_first_of(',');
    int const semicolon_index = string.find_first_of(';');
    if (((comma_index < semicolon_index) || (semicolon_index == -1)) && (comma_index > 0)) {
      vector.push_back(ConvertStringToNumber<TYPE>(TrimSpacesAndSquareBrackets(string.substr(0, comma_index))));
      string.erase(0, comma_index + 1);
    } else if (((semicolon_index < comma_index) || (comma_index == -1)) && (semicolon_index > 0)) {
      vector.push_back(ConvertStringToNumber<TYPE>(TrimSpacesAndSquareBrackets(string.substr(0, semicolon_index))));
      string.erase(0, semicolon_index + 1);
    } else {
      string.erase(0, 1);
    }
  }
  return vector;
}

template<typename TYPE>
std::string GetStringWithHighPrecision(const TYPE value) {
  std::ostringstream out;
  out.precision(16);
  out << value;
  return out.str();
}
