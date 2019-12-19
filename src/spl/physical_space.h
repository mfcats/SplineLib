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
  PhysicalSpace(const std::vector<spl::ControlPoint> &control_points,
                std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points);
  PhysicalSpace(const PhysicalSpace &physical_space);
  PhysicalSpace(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&other) noexcept = default;
  PhysicalSpace & operator=(const PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &rhs) = default;
  PhysicalSpace & operator=(PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &&rhs) noexcept = default;
  virtual ~PhysicalSpace() = default;

  bool AreEqual(const PhysicalSpace<PARAMETRIC_DIMENSIONALITY> &rhs,
                double tolerance = util::numeric_settings::GetEpsilon<double>()) const;

  virtual spl::ControlPoint GetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const;
  virtual spl::ControlPoint GetControlPoint(int index_1d) const;

  void SetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, const spl::ControlPoint &control_point,
                       int dimension = 0, int (*before)(int) = nullptr);
  void SetControlPoint(int index_1d, const spl::ControlPoint &control_point,
                       int dimension = 0, int (*before)(int) = nullptr);

  virtual void AddControlPoints(int number);
  virtual void RemoveControlPoints(int number);

  virtual int GetNumberOfControlPoints() const;
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetPointsPerDirection() const;
  std::array<int, PARAMETRIC_DIMENSIONALITY> GetMaximumPointIndexInEachDirection() const;

  void IncrementNumberOfPoints(int dimension);
  void DecrementNumberOfPoints(int dimension);

  void SetNumberOfPoints(int dimension, int value);

  virtual int GetDimension() const;

  virtual double GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> /*indices*/) const;
  virtual double GetWeight(int /*index_1d*/) const;

  double GetExpansion() const;

  double GetMaximumDistanceFromOrigin() const;

  std::vector<spl::ControlPoint> GetDividedControlPoints(int first, int length, int dimension);

 protected:
  int dimensionality_{};
  std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_{};
  std::vector<double> control_points_;
};

#include "src/spl/physical_space.inc"
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_PHYSICAL_SPACE_H_
