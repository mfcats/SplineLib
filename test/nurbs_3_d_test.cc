/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <numeric_settings.h>
#include "gmock/gmock.h"

#include "nurbs.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;
using testing::Return;
using ::testing::NiceMock;
using ::testing::Throw;
using ::testing::_;

class MockParameterSpace3d : public spl::ParameterSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 3>, std::array<ParamCoord, 3>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 3>, std::array<ParamCoord, 3>, std::array<int, 3>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 3>(std::array<ParamCoord, 3>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 3>));
};

class MockWeightedPhysicalSpace3d : public spl::WeightedPhysicalSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetWeight, double(std::array<int, 3>));
  MOCK_CONST_METHOD1(GetHomogenousControlPoint, baf::ControlPoint(std::array<int, 3>));
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 3>));
};

class MockPhysicalSpace3d : public spl::PhysicalSpace<3> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 3>));
};

void mock_weights3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
ON_CALL(*w_physical_space, GetWeight(_))
.WillByDefault(Return(1));
}

void mock_homogenous3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetHomogenousControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 0.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*w_physical_space, GetControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 2.0})));
}

void mock_weightedPhysicalSpace3d(const std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> &w_physical_space) {
  mock_weights3d(w_physical_space);
  mock_homogenous3d(w_physical_space);
}

void mock_physicalSpace3d(const std::shared_ptr<NiceMock<MockPhysicalSpace3d>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 0, 0}))
      .WillByDefault(Return(baf::ControlPoint({3.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 1, 0}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.5})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 0, 1}))
      .WillByDefault(Return(baf::ControlPoint({2.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{0, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({4.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 2.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 3>{2, 1, 1}))
      .WillByDefault(Return(baf::ControlPoint({5.0, 2.0})));
}

void set_get_basis_function_nurbs3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctions(_, std::array<ParamCoord, 3>{ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.5},
                                                                        ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.5},
                                                                        ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.5},
                                                                        ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.5},
                                                                        ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 0, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.045));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.09));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 0, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.045));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 1, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.36));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 1, 0},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.18));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 0, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.005));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 0, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.01));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 0, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.005));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{0, 1, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.02));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{1, 1, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.04));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 3>{2, 1, 1},
                                              std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                        ParamCoord{0.8},
                                                                        ParamCoord{0.1}}))
      .WillByDefault(Return(0.02));

}

void set_basis_function_derivative3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.5},
                                                                                  ParamCoord{0.5}},
                                                        std::array<int, 3>{1, 1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 0},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 0, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 0, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{0, 1, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 3>{2, 1, 1},
                                                        std::array<ParamCoord, 3>{ParamCoord{0.5},
                                                                                  ParamCoord{0.8},
                                                                                  ParamCoord{0.1}},
                                                        std::array<int, 3>{1, 1, 1}))
      .WillByDefault(Return(1.0));
}

void mock_parameterSpace_nurbs3d(const std::shared_ptr<NiceMock<MockParameterSpace3d>> &parameter_space) {
  set_get_basis_function_nurbs3d(parameter_space);
  set_basis_function_derivative3d(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault(Return(std::array<int, 3>{0, 0, 0}));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{1}));
  ON_CALL(*parameter_space, GetDegree(0))
      .WillByDefault(Return(Degree{2}));

}


class A3DNurbsWithAllWeights1 : public Test {
 public:
  A3DNurbsWithAllWeights1() :
      parameter_space_m(std::make_shared<NiceMock<MockParameterSpace3d>>()),
      w_physical_space_m(std::make_shared<NiceMock<MockWeightedPhysicalSpace3d>>()),
      physical_space_m(std::make_shared<NiceMock<MockPhysicalSpace3d>>()) {
    std::array<baf::KnotVector, 3> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}),
         baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}})};
    std::array<Degree, 3> degree = {Degree{2}, Degree{1}, Degree{1}};
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({3.0, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({5.0, 2.0}))
    };
    std::array<int, 3> number_of_points = {3, 2, 2};
    std::array<std::shared_ptr<baf::KnotVector>, 3> knot_vector_ptr =
        {std::make_shared<baf::KnotVector>(knot_vector[0]),
         std::make_shared<baf::KnotVector>(knot_vector[1]),
         std::make_shared<baf::KnotVector>(knot_vector[2])};

    parameter_space = std::make_shared<spl::ParameterSpace<3>>(knot_vector_ptr, degree);
    w_physical_space = std::make_shared<spl::WeightedPhysicalSpace<3>>(control_points, weights, number_of_points);
    physical_space = std::make_shared<spl::PhysicalSpace<3>>(control_points, number_of_points);
    /*
    spl::NURBSGenerator<3> nurbs_generator(w_physical_space, parameter_space);
    spl::BSplineGenerator<3> bspline_generator(physical_space, parameter_space);
    */
    spl::NURBSGenerator<3> nurbs_generator(w_physical_space_m, parameter_space_m);
    spl::BSplineGenerator<3> bspline_generator(physical_space_m, parameter_space_m);


    nurbs_ = std::make_unique<spl::NURBS<3>>(nurbs_generator);
    bspline_ = std::make_unique<spl::BSpline<3>>(bspline_generator);

  }

 protected:
  std::unique_ptr<spl::NURBS<3>> nurbs_;
  std::unique_ptr<spl::BSpline<3>> bspline_;
  std::shared_ptr<spl::ParameterSpace<3>> parameter_space;
  std::shared_ptr<spl::WeightedPhysicalSpace<3>> w_physical_space;
  std::shared_ptr<spl::PhysicalSpace<3>> physical_space;
  std::shared_ptr<NiceMock<MockParameterSpace3d>> parameter_space_m;
  std::shared_ptr<NiceMock<MockWeightedPhysicalSpace3d>> w_physical_space_m;
  std::shared_ptr<NiceMock<MockPhysicalSpace3d>> physical_space_m;
};

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_5And0_5AndDerivatives1And1And0) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space_m);
  mock_physicalSpace3d(physical_space_m);
  mock_weightedPhysicalSpace3d(w_physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 1, 0})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}, ParamCoord{0.5}},
                                                    {0},
                                                    {1, 1, 0})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And1And1) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space_m);
  mock_physicalSpace3d(physical_space_m);
  mock_weightedPhysicalSpace3d(w_physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}}, {0}, {1, 1, 1})[0],
              DoubleEq(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}},
                                                    {0},
                                                    {1, 1, 1})[0]));
}

TEST_F(A3DNurbsWithAllWeights1, ReturnsSameDerivativeAs3DBSplineFor0_5And0_8And0_1AndDerivatives1And2And1) { // NOLINT
  mock_parameterSpace_nurbs3d(parameter_space_m);
  mock_physicalSpace3d(physical_space_m);
  mock_weightedPhysicalSpace3d(w_physical_space_m);
  ASSERT_THAT(nurbs_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}}, {0}, {1, 2, 1})[0],
              DoubleNear(bspline_->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.8}, ParamCoord{0.1}},
                                                      {0},
                                                      {1, 2, 1})[0],
                         util::NumericSettings<double>::kEpsilon()));
}
