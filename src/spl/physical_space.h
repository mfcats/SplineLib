/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_PHYSICAL_SPACE_H_
#define SRC_SPL_PHYSICAL_SPACE_H_

#include <stdexcept>
#include <vector>

#include "src/spl/control_point.h"
#include "src/util/multi_index_handler.h"
#include "src/util/numeric_settings.h"
#include "src/util/vector_utils.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class PhysicalSpace {
 public:
  PhysicalSpace() = default;
  PhysicalSpace(std::vector<ControlPoint> const &control_points,
                std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points);
  PhysicalSpace(PhysicalSpace const &physical_space) = default;
  PhysicalSpace(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&other) noexcept = default;
  PhysicalSpace & operator=(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> const &rhs) = default;
  PhysicalSpace & operator=(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&rhs) noexcept = default;
  virtual ~PhysicalSpace() = default;

  virtual int GetDimensionality() const;
  int GetNumberOfPointsForDimension(Dimension const &dimension) const;
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetNumberOfPointsPerDirection() const;
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetMaximumPointIndexPerDirection() const;
  virtual int GetTotalNumberOfControlPoints() const;

  virtual ControlPoint GetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const;
  virtual ControlPoint GetControlPoint(int index_1d) const;
  std::vector<double> GetControlPoints() const;

  virtual Weight GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> const &/*indices*/) const;
  virtual Weight GetWeight(int /*index_1d*/) const;
  std::vector<Weight> GetWeights() const {
    return std::vector<Weight>(total_number_of_points, Weight{1.0});
  }

  void SetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, ControlPoint const &control_point,
                       Dimension const &dimension = Dimension{0}, int (*before)(int) = nullptr);
  void SetControlPoint(int index_1d, ControlPoint const &control_point,
                       Dimension const &dimension = Dimension{0}, int (*before)(int) = nullptr);

  double GetExpansion() const;
  double GetMaximumDistanceFromOrigin() const;

  std::vector<ControlPoint> SplitControlPoints(Dimension const &dimension, int offset, int length);

  bool AreEqual(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> const &rhs,
                Tolerance const &tolerance = Tolerance{util::numeric_settings::GetEpsilon<double>()}) const;

  // TODO(all): The following methods should become protected.
  void IncrementNumberOfPoints(int dimension);
  void DecrementNumberOfPoints(int dimension);

  virtual void AddControlPoints(int number);
  virtual void RemoveControlPoints(int number);

  void SetNumberOfPoints(int dimension, int value);

 protected:
  int dimensionality_{};
  std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_per_dimension_{};
  int total_number_of_points;
  std::vector<double> control_points_;
};

#include "src/spl/physical_space.inc"
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_PHYSICAL_SPACE_H_
