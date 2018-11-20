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
#include <memory>
#include <sstream>
#include <vector>

#include "alias.h"
#include "basis_function.h"
#include "basis_function_factory.h"
#include "knot_vector.h"

namespace spl {
template<int DIM>
class ParameterSpace {
 public:
  ParameterSpace() = default;

  ParameterSpace(const std::array<std::shared_ptr<baf::KnotVector>, DIM> &knot_vector, std::array<Degree, DIM> degree) :
      knot_vector_(knot_vector), degree_(degree) {
    ThrowIfKnotVectorDoesNotStartAndEndWith();
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

  std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator GetFirstNonZeroKnot(int direction,
                                                                                       ParamCoord param_coord) const {
    return basis_functions_[direction].begin() + knot_vector_[direction]->GetKnotSpan(param_coord).get()
        - degree_[direction].get();
  }

  virtual Degree GetDegree(int direction) const {
    return degree_[direction];
  }

  std::array<Degree, DIM> GetDegrees() const {
    return degree_;
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int direction) const {
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

  std::array<std::vector<ParamCoord>, DIM> GetKnots() const {
    std::array<std::vector<ParamCoord>, DIM> knots;
    for (int i = 0; i < DIM; ++i) {
      std::vector<ParamCoord> temp = knot_vector_[i]->GetKnots();
      for (auto &j : temp) {
        knots[i].emplace_back(j);
      }
    }
    return knots;
  }

  void InsertKnot(ParamCoord knot, int dimension) {
    knot_vector_[dimension]->InsertKnot(knot);
    for (int i = 0; i < DIM; i++) {
      basis_functions_[i].erase(basis_functions_[i].begin(), basis_functions_[i].end());
      basis_functions_[i].reserve(knot_vector_[i]->GetNumberOfKnots() - degree_[i].get() - 1);
      for (int j = 0; j < (static_cast<int>(knot_vector_[i]->GetNumberOfKnots()) - degree_[i].get() - 1); ++j) {
        basis_functions_[i].emplace_back(baf::BasisFunctionFactory::CreateDynamic(*(knot_vector_[i]),
                                                                                  KnotSpan{j},
                                                                                  degree_[i]));
      }
    }
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

  KnotVectors<DIM> knot_vector_;
  std::array<Degree, DIM> degree_;
  std::array<std::vector<std::shared_ptr<baf::BasisFunction>>, DIM> basis_functions_;
};
}  // namespace spl

#endif  // SRC_SPL_PARAMETER_SPACE_H_
