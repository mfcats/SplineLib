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
  MultiIndexHandler(std::array<int, DIM> &multi_index_length)
      : multi_index_length(multi_index_length), current_multi_index_value({0}) {}

  int operator[](int i) {
    return current_multi_index_value[i];
  }

  MultiIndexHandler &operator++() {
    for (int i = 0; i < DIM; ++i) {
      if (current_multi_index_value[i] == multi_index_length[i] - 1) {
        current_multi_index_value[i] = 0;
      } else {
        current_multi_index_value[i]++;
        break;
      }
    }
    return *this;
  }

  MultiIndexHandler operator++(int) {
    MultiIndexHandler result(*this);
    ++(*this);
    return result;
  }

  void SetIndices(std::array<int, DIM> &indices){
	  for (int i = 0; i < DIM; ++i) {
        current_multi_index_value[i] = indices[i];
	  }
  }

  int Get1DIndex(){
    int index_1d = 0;
    int temp;
	  for (int i = 0; i < DIM; ++i){
        temp = current_multi_index_value[i];
		  for(int j = i-1; j >= 0; --j){
            temp *= multi_index_length[j];
		  }
        index_1d += temp;
	  }
    return index_1d;
  }

 private:
  std::array<int, DIM> multi_index_length;
  std::array<int, DIM> current_multi_index_value;
};

#endif //SPLINELIB_MULTI_INDEX_HANDLER_H
