/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_IGA_MATRIX_H_
#define SRC_IGA_MATRIX_H_

#include <array>

namespace iga {
class Matrix {
 public:
  explicit Matrix(int rowCount, int colCount) : rowCount_(rowCount), colCount_(colCount) {
    matrix_ = new double[rowCount_ * colCount_]();
  }

  ~Matrix() {
    delete[] matrix_;
  }

  void WriteToMatrix(int row, int col, double number) {
    matrix_[row * colCount_ + col] = number;
  }

  void AddToMatrixEntry(int row, int col, double number) {
    matrix_[row * colCount_ + col] += number;
  }

  double GetMatrixEntry(int row, int col) {
    return matrix_[row * colCount_ + col];
  }

  int GetRowCount() {
    return rowCount_;
  }

  int GetColCount() {
    return colCount_;
  }

 private:
  int rowCount_;
  int colCount_;
  double *matrix_;
};
}  // namespace iga

#endif  // SRC_IGA_MATRIX_H_
