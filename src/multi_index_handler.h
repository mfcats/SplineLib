/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_MULTI_INDEX_HANDLER_H
#define SPLINELIB_MULTI_INDEX_HANDLER_H

#include "multi_index_handler.h"
#include <array>

template<int DIM>
class MultiIndexHandler {
 public:
  MultiIndexHandler(std::array<int, DIM> &lastKnotOffset)
      : lastKnotOffset(lastKnotOffset), currentMultiIndex({0}) {}

  int operator[](int i) {
    return currentMultiIndex[i];
  }

  MultiIndexHandler &operator++() {
    for (int i = 0; i < DIM; ++i) {
      if (currentMultiIndex[i] == lastKnotOffset[i]) {
        currentMultiIndex[i] = 0;
      } else {
        currentMultiIndex[i]++;
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
  std::array<int, DIM> lastKnotOffset;
  std::array<int, DIM> currentMultiIndex;
};

#endif //SPLINELIB_MULTI_INDEX_HANDLER_H
