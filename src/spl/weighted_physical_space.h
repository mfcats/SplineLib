/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
#define SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_

#include <utility>
#include <vector>

#include "src/spl/physical_space.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class WeightedPhysicalSpace : public PhysicalSpace<PARAMETRIC_DIMENSIONALITY> {
 public:
  WeightedPhysicalSpace() = default;
  WeightedPhysicalSpace(std::vector<ControlPoint> const &control_points, std::vector<Weight> const &weights,
                        std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points);
  WeightedPhysicalSpace(WeightedPhysicalSpace const &other) = default;
  WeightedPhysicalSpace(WeightedPhysicalSpace &&other) noexcept = default;
  WeightedPhysicalSpace & operator=(WeightedPhysicalSpace const &rhs) = default;
  WeightedPhysicalSpace & operator=(WeightedPhysicalSpace &&rhs) noexcept = default;
  ~WeightedPhysicalSpace() override = default;

  Weight GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const override;
  Weight GetWeight(int index_1d) const override;
  std::vector<Weight> GetWeights() const override;

  void SetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, Weight const &weight,
                 Dimension const &dimension = Dimension{0}, int (*before)(int) = nullptr);
  void SetWeight(int index_1d, Weight const &weight, Dimension const &dimension = Dimension{0},
                 int (*before)(int) = nullptr);

  virtual ControlPoint GetHomogeneousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const;
  virtual ControlPoint GetHomogeneousControlPoint(int index_1d) const;

  void SetHomogeneousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
                                  ControlPoint const &homogeneous_control_point);
  void SetWeightedControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
                               ControlPoint const &control_point, Weight const &weight);

  double GetMinimumWeight() const;

  std::vector<Weight> SplitWeights(Dimension dimension, int offset, int length);

  bool AreEqual(WeightedPhysicalSpace const &rhs,
                Tolerance const &tolerance = Tolerance{util::numeric_settings::GetEpsilon<double>()}) const;

  // TODO(all): The arguments should be a dimension and an index for this dimension or the method should be private.
  void AddControlPoints(int number) override;
  void DoubleControlPointSlice(Dimension const &dimension, int index) override;
  void RemoveControlPoints(int number) override;

 private:
  std::vector<Weight> weights_;
};

#include "src/spl/weighted_physical_space.inc"
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
