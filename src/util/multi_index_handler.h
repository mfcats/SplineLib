/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_MULTI_INDEX_HANDLER_H_
#define SRC_UTIL_MULTI_INDEX_HANDLER_H_

#include <array>
#include <utility>

#include "named_type.h"
#include "stl_container_access.h"

namespace splinelib::src::util {
// MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>s map a PARAMETRIC_DIMENSIONALITY-dimensional tensor-product array to an
// one-dimensional array by mapping a PARAMETRIC_DIMENSIONALITY-dimensional multi-index to an integer (and vice versa).
// It also provides "iterator-like" functionality.
// Example ((mixed) derivatives up to fourth order with respect to three variables (PARAMETRIC_DIMENSIONALITY = 3)):
//   array<int, 3> maximum_order_derivatives = {3, 2, 2};
//   MultiIndexHandler<2> derivative_handler(maximum_order_derivatives);  // Initial value of multi-index: {0, 0, 0}.
//   ++derivative_handler;  // Move from {0, 0, 0} to {1, 0, 0}.
//   array<int, 3> multi_index_d3/dx2dy = {2, 1, 0};
//   derivative_handler.SetCurrentIndex(multi_index_d3/dx2dy);  // Set multi-index to {2, 1, 0}.
//   1d_index_d3/dx2dy = Get1DIndex();  // Store 5 in 1d_index_d3/dx2dy.
//   ++derivative_handler;  // Move from {2, 1, 0} to {0, 0, 1}.
template<int PARAMETRIC_DIMENSIONALITY>
class MultiIndexHandler {
 public:
  MultiIndexHandler() = default;
  explicit MultiIndexHandler(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index_length);
  MultiIndexHandler(MultiIndexHandler const &other) = default;
  MultiIndexHandler(MultiIndexHandler &&other) noexcept;
  MultiIndexHandler & operator=(MultiIndexHandler const &rhs) = default;
  MultiIndexHandler & operator=(MultiIndexHandler &&rhs) noexcept;
  ~MultiIndexHandler() = default;

  int operator[](Dimension const &dimension) const;
  MultiIndexHandler & operator++();
  MultiIndexHandler operator++(int);
  MultiIndexHandler & operator--();
  MultiIndexHandler operator--(int);
  // Jumps count multi-index values forwards.
  MultiIndexHandler & operator+(int count);
  // Jumps count multi-index values backwards.
  MultiIndexHandler & operator-(int count);

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetCurrentIndex() const;
  int GetCurrent1DIndex() const;
  static int Get1DIndex(std::array<int, PARAMETRIC_DIMENSIONALITY> const &length,
                        std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index);
  int Get1DIndex(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const;
  void SetCurrentIndex(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index);
  void SetCurrentIndex(int index_1d);

  int Get1DLength() const;
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetComplementaryIndex() const;

  int CollapseDimension(Dimension const &dimension) const;

 private:
  static int GetCurrentSliceSize(std::array<int, PARAMETRIC_DIMENSIONALITY> const &length,
                                 Dimension const &maximum_dimension);

  std::array<int, PARAMETRIC_DIMENSIONALITY> multi_index_length_;
  std::array<int, PARAMETRIC_DIMENSIONALITY> current_multi_index_value_;
};

#include "multi_index_handler.inc"
}  // namespace splinelib::src::util

#endif  // SRC_UTIL_MULTI_INDEX_HANDLER_H_
