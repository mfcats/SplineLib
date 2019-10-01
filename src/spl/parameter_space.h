/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PARAMETER_SPACE_H_
#define SRC_SPL_PARAMETER_SPACE_H_

#include <array>
#include <algorithm>
#include <memory>
#include <sstream>
#include <vector>

#include "alias.h"
#include "basis_function.h"
#include "basis_function_factory.h"
#include "knot_vector.h"
#include "numeric_settings.h"

namespace spl {
template<int DIM>
class ParameterSpace {
 public:
  ParameterSpace() = default;

  ParameterSpace(const KnotVectors<DIM> &knot_vector, std::array<Degree, DIM> degree)
      : knot_vector_(knot_vector), degree_(degree) {
    ThrowIfKnotVectorDoesNotStartAndEndWith();
    baf::BasisFunctionFactory factory;
    for (int i = 0; i < DIM; i++) {
      basis_functions_[i].reserve(knot_vector_[i]->GetNumberOfKnots() - degree_[i].get() - 1);
      for (int j = 0; j < (static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) - degree_[i].get() - 1); ++j) {
        basis_functions_[i].emplace_back(factory.CreateDynamic(*(knot_vector_[i]), KnotSpan{j}, degree_[i]));
      }
    }
  }

  ParameterSpace(const ParameterSpace<DIM> &parameter_space) {
    degree_ = parameter_space.degree_;
    for (int i = 0; i < DIM; ++i) {
      baf::KnotVector knot_vector = *(parameter_space.GetKnotVector(i));
      knot_vector_[i] = std::make_shared<baf::KnotVector>(knot_vector);
    }
    baf::BasisFunctionFactory factory;
    for (int i = 0; i < DIM; i++) {
      basis_functions_[i].reserve(knot_vector_[i]->GetNumberOfKnots() - degree_[i].get() - 1);
      for (int j = 0; j < (static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) - degree_[i].get() - 1); ++j) {
        basis_functions_[i].emplace_back(factory.CreateDynamic(*(knot_vector_[i]), KnotSpan{j}, degree_[i]));
      }
    }
  }

  virtual ~ParameterSpace() = default;

  std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParamCoord param_coord) const {
    auto first_non_zero = GetFirstNonZeroKnot(direction, param_coord);
    std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_[direction].get()) + 1, 0.0);
    for (int i = 0; i < degree_[direction].get() + 1; ++i) {
      basis_function_values[i] = (*first_non_zero)->Evaluate(param_coord);
      ++first_non_zero;
    }
    return basis_function_values;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
                                                                 ParamCoord param_coord,
                                                                 int derivative) const {
    auto first_non_zero = GetFirstNonZeroKnot(direction, param_coord);
    std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_[direction].get()) + 1, 0.0);
    for (int i = 0; i < degree_[direction].get() + 1; ++i) {
      basis_function_values[i] = (*first_non_zero)->EvaluateDerivative(param_coord, Derivative{derivative});
      ++first_non_zero;
    }
    return basis_function_values;
  }

  bool AreEqual(const ParameterSpace<DIM> &rhs, double tolerance = util::NumericSettings<double>::kEpsilon()) const {
    return std::equal(degree_.begin(), degree_.end(), rhs.degree_.begin(), rhs.degree_.end(),
                      [&](Degree degree_a, Degree degree_b) {
                        return util::NumericSettings<double>::AreEqual(degree_a.get(), degree_b.get());
                      })
        && std::equal(knot_vector_.begin(), knot_vector_.end(), rhs.knot_vector_.begin(), rhs.knot_vector_.end(),
                      [&](std::shared_ptr<baf::KnotVector> kv_a, std::shared_ptr<baf::KnotVector> kv_b) {
                        return kv_a->AreEqual(*kv_b.get(), tolerance);
                      });
  }

  std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator
  GetFirstNonZeroKnot(int direction, ParamCoord param_coord) const {
    return basis_functions_[direction].begin() + knot_vector_[direction]->GetKnotSpan(param_coord).get()
        - degree_[direction].get();
  }

  virtual Degree GetDegree(int direction) const {
    return degree_[direction];
  }

  virtual std::shared_ptr<baf::KnotVector> GetKnotVector(int direction) const {
    return knot_vector_[direction];
  }

  virtual double GetBasisFunctions(std::array<int, DIM> indices, std::array<ParamCoord, DIM> param_coord) const {
    double value = 1;
    for (int i = 0; i < DIM; ++i) {
      value *= basis_functions_[i][indices[i]]->Evaluate(param_coord[i]);
    }
    return value;
  }

  virtual double GetBasisFunctionDerivatives(std::array<int, DIM> indices,
                                             std::array<ParamCoord, DIM> param_coord,
                                             std::array<int, DIM> derivative) const {
    double value = 1;
    for (int i = 0; i < DIM; ++i) {
      value *= basis_functions_[i][indices[i]]->EvaluateDerivative(param_coord[i], Derivative{derivative[i]});
    }
    return value;
  }

  virtual std::array<int, DIM> GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const {
    std::array<int, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = GetKnotVector(i)->GetKnotSpan(param_coord[i]).get() - GetDegree(i).get();
    }
    return first_non_zero;
  }

  virtual void ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, DIM> param_coord) const {
    for (int dim = 0; dim < DIM; dim++) {
      if (!this->GetKnotVector(dim)->IsInKnotVectorRange(param_coord[dim])) {
        std::stringstream message;
        message << "The parametric coordinate " << param_coord[dim].get() << " is outside the knot vector range from "
                << GetKnotVector(dim)->GetKnot(0).get() << " to " << GetKnotVector(dim)->GetLastKnot().get() << ".";
        throw std::range_error(message.str());
      }
    }
  }

  double GetKnotVectorRange(int direction) const {
    return GetKnotVector(direction)->GetLastKnot().get() - GetKnotVector(direction)->GetKnot(0).get();
  }

  void InsertKnot(ParamCoord knot, int dimension) {
    knot_vector_[dimension]->InsertKnot(knot);
    RecreateBasisFunctions();
  }

  void RemoveKnot(ParamCoord knot, int dimension) {
    knot_vector_[dimension]->RemoveKnot(knot);
    RecreateBasisFunctions();
  }

  void ElevateDegree(int dimension) {
    degree_[dimension] = Degree{degree_[dimension].get() + 1};
    RecreateBasisFunctions();
  }

  void ReduceDegree(int dimension) {
    degree_[dimension] = Degree{degree_[dimension].get() - 1};
    RecreateBasisFunctions();
  }

  void IncrementMultiplicityOfAllKnots(int dim) {
    std::vector<ParamCoord> unique_knots(GetKnotVector(dim)->GetNumberOfDifferentKnots());
    auto last = std::unique_copy(GetKnotVector(dim)->begin(), GetKnotVector(dim)->end(), unique_knots.begin());
    unique_knots.resize(std::distance(unique_knots.begin(), last));
    for (auto &knot : unique_knots) {
      InsertKnot(knot, dim);
    }
  }

  void DecrementMultiplicityOfAllKnots(int dim) {
    std::vector<ParamCoord> unique_knots(GetKnotVector(dim)->GetNumberOfDifferentKnots());
    auto last = std::unique_copy(GetKnotVector(dim)->begin(), GetKnotVector(dim)->end(), unique_knots.begin());
    unique_knots.resize(std::distance(unique_knots.begin(), last));
    for (auto &knot : unique_knots) {
      RemoveKnot(knot, dim);
    }
  }

  std::array<KnotVectors<DIM>, 2> GetDividedKnotVectors(ParamCoord param_coord, int dimension) const {
    auto knot_span = GetKnotVector(dimension)->GetKnotSpan(param_coord).get();
    std::array<int, 2> first_knot = {0, knot_span - GetDegree(dimension).get()};
    std::array<int, 2> last_knot = {knot_span + 1, static_cast<int>(knot_vector_[dimension]->GetNumberOfKnots())};
    std::array<KnotVectors<DIM>, 2> new_knot_vectors;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < DIM; ++j) {
        new_knot_vectors[i][j] = GetKnotVector(j);
      }
      std::array<std::vector<ParamCoord>, 2> new_knots;
      for (int j = first_knot[i]; j < last_knot[i]; ++j) {
        new_knots[i].push_back(knot_vector_[dimension]->GetKnot(j));
      }
      new_knot_vectors[i][dimension] = std::make_shared<baf::KnotVector>(baf::KnotVector(new_knots[i]));
    }
    return new_knot_vectors;
  }

 private:
  void ThrowIfKnotVectorDoesNotStartAndEndWith() {
    for (int i = 0; i < DIM; i++) {
      if (static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) < 2 * degree_[i].get() + 2) {
        throw std::runtime_error("There have to be at least 2p + 2 knots.");
      }
      for (int j = 1; j < degree_[i].get() + 1; j++) {
        if (knot_vector_[i]->GetKnot(0).get() != knot_vector_[i]->GetKnot(j).get()) {
          throw std::runtime_error("The first knot must have multiplicity p+1.");
        }
      }
      for (int j = static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) - 2;
           j > static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) - degree_[i].get() - 2; j--) {
        if (knot_vector_[i]->GetLastKnot().get() != knot_vector_[i]->GetKnot(j).get()) {
          throw std::runtime_error("The last knot must have multiplicity p+1.");
        }
      }
    }
  }

  void RecreateBasisFunctions() {
    for (int current_dim = 0; current_dim < DIM; ++current_dim) {
      auto &basis_functions_current_dim = basis_functions_[current_dim];
      auto &knot_vector_current_dim = knot_vector_[current_dim];
      auto &degree_current_dim = degree_[current_dim];
      basis_functions_current_dim.erase(basis_functions_current_dim.begin(), basis_functions_current_dim.end());
      basis_functions_current_dim.reserve(knot_vector_current_dim->GetNumberOfKnots() - degree_current_dim.get() - 1);
      for (int current_basis_function = 0;
           current_basis_function <
               (static_cast<int>(knot_vector_current_dim->GetNumberOfKnots()) - degree_current_dim.get() - 1);
           ++current_basis_function) {
        basis_functions_current_dim.emplace_back(baf::BasisFunctionFactory::CreateDynamic(
            *knot_vector_current_dim, KnotSpan{current_basis_function}, degree_current_dim));
      }
    }
  }

  KnotVectors<DIM> knot_vector_;
  std::array<Degree, DIM> degree_;
  std::array<std::vector<std::shared_ptr<baf::BasisFunction>>, DIM> basis_functions_;
};
}  // namespace spl

#endif  // SRC_SPL_PARAMETER_SPACE_H_
