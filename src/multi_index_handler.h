/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_MULTI_INDEX_HANDLER_H_
#define SRC_MULTI_INDEX_HANDLER_H_

#include <array>

template<int DIM>
class MultiIndexHandler {
 public:
  explicit MultiIndexHandler(const std::array<int, DIM> &maximum_multi_index_value) : maximum_multi_index_value_(
      maximum_multi_index_value), current_multi_index_value_({0}) {}

  int operator[](int i) {
    return current_multi_index_value_[i];
  }

  MultiIndexHandler &operator++() {
    for (int i = 0; i < DIM; ++i) {
      if (current_multi_index_value_[i] == maximum_multi_index_value_[i]) {
        current_multi_index_value_[i] = 0;
      } else {
        current_multi_index_value_[i]++;
        break;
      }
    }
    return *this;
  }

  const MultiIndexHandler operator++(int) {
    MultiIndexHandler result(*this);
    ++(*this);
    return result;
  }

 private:
  std::array<int, DIM> maximum_multi_index_value_;
  std::array<int, DIM> current_multi_index_value_;
};

#endif  // SRC_MULTI_INDEX_HANDLER_H_
