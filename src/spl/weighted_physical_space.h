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
  WeightedPhysicalSpace(std::vector<spl::ControlPoint> const &control_points, std::vector<double> const &weights,
                        std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points);
  WeightedPhysicalSpace(WeightedPhysicalSpace const &other);
  WeightedPhysicalSpace(WeightedPhysicalSpace &&other) noexcept = default;
  WeightedPhysicalSpace & operator=(WeightedPhysicalSpace const &rhs) = default;
  WeightedPhysicalSpace & operator=(WeightedPhysicalSpace &&rhs) noexcept = default;
  ~WeightedPhysicalSpace() override = default;

  bool AreEqual(WeightedPhysicalSpace const &rhs, double tolerance = util::numeric_settings::GetEpsilon<double>()) const;

  virtual ControlPoint GetHomogenousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const;

  virtual ControlPoint GetHomogenousControlPoint(int index_1d) const;

  Weight GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const override;

  Weight GetWeight(int index_1d) const override;

  double GetMinimumWeight() const;

  void SetWeightedControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
                               ControlPoint const &control_point, double weight);

  void SetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, double weight, int dimension = 0,
                 int (*before)(int) = nullptr);

  void SetWeight(int index_1d, double weight, int dimension = 0,
                 int (*before)(int) = nullptr);

  void AddControlPoints(int number) final;

  void RemoveControlPoints(int number) final;

  virtual std::vector<double> GetWeights() const;  // final;  Weight instead of double

  std::vector<double> GetDividedWeights(int first, int length, int dimension);

 private:
  std::vector<double> weights_;
};

#include "src/spl/weighted_physical_space.inc"
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_WEIGHTED_PHYSICAL_SPACE_H_
