/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_MATRIX_UTILS_H_
#define SRC_IGA_MATRIX_UTILS_H_

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace iga {
class MatrixUtils {
 public:
  template<typename T>
  static void PrintMatrix(std::vector<std::vector<T>> matrix) {
    std::cout << std::endl;
    for (int i = 0; i < matrix.size(); ++i) {
      for (int j = 0; j < matrix[0].size(); ++j) {
        uint64_t max_length = std::to_string(matrix[matrix.size() - 1][matrix[0].size() - 1]).length();
        std::string temp = std::to_string(matrix[i][j]);
        for (uint64_t i = 0; i < max_length - temp.length(); ++i) {
          temp.insert(0, " ");
        }
        std::cout << temp << "   ";
      }
      std::cout << std::endl;
    }
  }

  static std::array<std::array<double, 2>, 2> Get2By2Inverse(std::array<std::array<double, 2>, 2> matrix) {
    std::array<std::array<double, 2>, 2> inverse;
    double factor = 1 / Get2By2Determinant(matrix);
    inverse[0][0] = factor * matrix[1][1];
    inverse[0][1] = factor * matrix[0][1] * (-1);
    inverse[1][0] = factor * matrix[1][0] * (-1);
    inverse[1][1] = factor * matrix[0][0];
    return inverse;
  }

  static double Get2By2Determinant(std::array<std::array<double, 2>, 2> matrix) {
    return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
  }
};
}  // namespace iga

#endif  // SRC_IGA_MATRIX_UTILS_H_
