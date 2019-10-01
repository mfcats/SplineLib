/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_MULTI_INDEX_HANDLER_H_
#define SRC_UTIL_MULTI_INDEX_HANDLER_H_

#include <array>
#include <vector>

namespace util {
template<int DIM>
class MultiIndexHandler {
 public:
  MultiIndexHandler() = default;
  explicit MultiIndexHandler(const std::array<int, DIM> &multi_index_length) : multi_index_length_(multi_index_length) {
    for (auto &index : current_multi_index_value_) {
      index = 0;
    }
  }

  int operator[](int i) {
#ifdef DEBUG
    return current_multi_index_value_.at(i);
#else
    return current_multi_index_value_[i];
#endif
  }

  MultiIndexHandler &operator++() {
    for (int i = 0; i < DIM; ++i) {
      if (current_multi_index_value_[i] == multi_index_length_[i] - 1 || multi_index_length_[i] == 0) {
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

  MultiIndexHandler &operator--() {
    for (int i = 0; i < DIM; ++i) {
      if (current_multi_index_value_[i] == 0) {
        current_multi_index_value_[i] = multi_index_length_[i] - 1;
      } else {
        current_multi_index_value_[i]--;
        break;
      }
    }
    return *this;
  }

  MultiIndexHandler &operator+(int i) {
    for (int j = 0; j < i; ++j) {
      ++(*this);
    }
    return *this;
  }

  MultiIndexHandler &operator-(int i) {
    for (int j = 0; j < i; ++j) {
      --(*this);
    }
    return *this;
  }

  void SetIndices(const std::array<int, DIM> &indices) {
    for (int i = 0; i < DIM; ++i) {
      current_multi_index_value_[i] = indices[i];
    }
  }

  void Set1DIndex(int index) {
    bool reached = false;
    while (!reached) {
      (*this)++;
      if (this->Get1DIndex() == index) reached = true;
    }
  }

  std::array<int, DIM> GetIndices() const {
    return current_multi_index_value_;
  }

  std::array<int, DIM> GetDifferenceIndices() const {
    std::array<int, DIM> indices;
    for (int i = 0; i < DIM; ++i) {
      indices[i] = multi_index_length_[i] - current_multi_index_value_[i] - 1;
    }
    return indices;
  }

  int Get1DIndex() const {
    return Get1DIndex(multi_index_length_, current_multi_index_value_);
  }

  int Get1DIndex(const std::array<int, DIM> &indices) const {
    return Get1DIndex(multi_index_length_, indices);
  }

  int Get1DIndex(const std::array<int, DIM> &length, const std::array<int, DIM> &indices) const {
    int index_1d = 0;
    for (int i = 0; i < DIM; ++i) {
      int temp = indices[i];
      for (int j = i - 1; j >= 0; --j) {
        temp *= length[j];
      }
      index_1d += temp;
    }
    return index_1d;
  }

  int Get1DLength() const {
    int length = 1;
    for (int i = 0; i < DIM; ++i) {
      length *= multi_index_length_[i];
    }
    return length;
  }

  int ExtractDimension(int dimension) const {
    std::array<int, DIM> indices, length;
    for (int m = 0; m < DIM; ++m) {
      if (m < dimension) {
        indices[m] = current_multi_index_value_[m];
        length[m] = multi_index_length_[m];
      } else if (m == DIM - 1) {
        indices[m] = 0;
        length[m] = 0;
      } else if (m >= dimension) {
        indices[m] = current_multi_index_value_[m + 1];
        length[m] = multi_index_length_[m + 1];
      }
    }
    return Get1DIndex(length, indices);
  }

 private:
  std::array<int, DIM> multi_index_length_;
  std::array<int, DIM> current_multi_index_value_;
};
}  // namespace util

#endif  // SRC_UTIL_MULTI_INDEX_HANDLER_H_
