/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <array>
#include <numeric>

#include "gmock/gmock.h"

#include "b_spline.h"
#include "one_point_gauss_legendre.h"
#include "two_point_gauss_legendre.h"
#include "three_point_gauss_legendre.h"
#include "four_point_gauss_legendre.h"
#include "five_point_gauss_legendre.h"
#include "numeric_settings.h"
#include "b_spline_generator.h"

using testing::Test;
using testing::DoubleEq;
using ::testing::Return;
using ::testing::Throw;
using ::testing::NiceMock;
using ::testing::_;

using namespace testing;

class Mock2dParameterSpace : public spl::ParameterSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetDegree, Degree(int));
  MOCK_CONST_METHOD2(GetBasisFunctions, double(std::array<int, 2>, std::array<ParamCoord, 2>));
  MOCK_CONST_METHOD3(GetBasisFunctionDerivatives,
                     double(std::array<int, 2>, std::array<ParamCoord, 2>, std::array<int, 2>));
  MOCK_CONST_METHOD1(GetArrayOfFirstNonZeroBasisFunctions, std::array<int, 2>(std::array<ParamCoord, 2>));
  MOCK_CONST_METHOD1(ThrowIfParametricCoordinateOutsideKnotVectorRange, void(std::array<ParamCoord, 2>));
};

class Mock2dPhysicalSpace : public spl::PhysicalSpace<2> {
 public:
  MOCK_CONST_METHOD1(GetControlPoint, baf::ControlPoint(std::array<int, 2>));
};

void mock_2dphysicalSpace(const std::shared_ptr<NiceMock<Mock2dPhysicalSpace>> &physical_space) {
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 0}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 0}))
      .WillByDefault(Return(baf::ControlPoint({0.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 0}))
      .WillByDefault(Return(baf::ControlPoint({1.0, -1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 1}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 1}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 0.0, 1.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 1}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 0.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{0, 2}))
      .WillByDefault(Return(baf::ControlPoint({-1.0, 1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{1, 2}))
      .WillByDefault(Return(baf::ControlPoint({0.0, 1.0, 0.0})));
  ON_CALL(*physical_space, GetControlPoint(std::array<int, 2>{2, 2}))
      .WillByDefault(Return(baf::ControlPoint({1.0, 1.0, 0.0})));
}

void set_get_basis_function(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {

  // Using wildcards
  ON_CALL(*parameter_space, GetBasisFunctions(_, _))
      .WillByDefault(Return(0.0));

  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.0}}))
      .WillByDefault(Return(1.0));

  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1}, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}}))
      .WillByDefault(Return(0.111109));

  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.33333}, ParamCoord{0.0}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0}, std::array<ParamCoord, 2>{ParamCoord{0.33333}, ParamCoord{0.0}}))
      .WillByDefault(Return(0.444449));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0}, std::array<ParamCoord, 2>{ParamCoord{0.33333}, ParamCoord{0.0}}))
      .WillByDefault(Return(0.111109));

  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.125));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2}, std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}}))
      .WillByDefault(Return(0.0625));

  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 0}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.0351562));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 0}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 0}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.316406));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 1}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 1}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.140625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 1}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.210938));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{0, 2}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.00390625));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{1, 2}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.0234375));
  ON_CALL(*parameter_space, GetBasisFunctions(std::array<int, 2>{2, 2}, std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}}))
      .WillByDefault(Return(0.0351562));

}

void set_basis_function_derivative1(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(_, _, _)).WillByDefault(Return(0.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                      std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.0}},
                                      std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.0}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.0}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(2.0));

  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.888898));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.888898));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.888884));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.888884));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.222218));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.222218));

  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-1.33334));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.66668));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.0}, ParamCoord{0.33333}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.66666));

  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.5));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.25));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.5}, ParamCoord{0.5}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.25));

  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.28125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.84375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.03125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(-0.0625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 0}))
      .WillByDefault(Return(0.09375));

  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.09375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(-0.84375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.0625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.375));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.5625));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.03125));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.1875));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{0, 1}))
      .WillByDefault(Return(0.28125));
}

void set_basis_function_derivative2(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 0},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(3.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(4.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 1},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-6.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{0, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-1.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{1, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(-2.0));
  ON_CALL(*parameter_space, GetBasisFunctionDerivatives(std::array<int, 2>{2, 2},
                                                        std::array<ParamCoord, 2>{ParamCoord{0.75}, ParamCoord{0.25}},
                                                        std::array<int, 2>{1, 2}))
      .WillByDefault(Return(3.0));
}

void mock_2dparameterSpace(const std::shared_ptr<NiceMock<Mock2dParameterSpace>> &parameter_space) {
  set_basis_function_derivative1(parameter_space);
  set_basis_function_derivative2(parameter_space);
  set_get_basis_function(parameter_space);
  ON_CALL(*parameter_space, GetArrayOfFirstNonZeroBasisFunctions(_))
      .WillByDefault((Return(std::array<int, 2>{0, 0})));
  ON_CALL(*parameter_space, GetDegree(_))
      .WillByDefault(Return(Degree{2}));
}

/* 2-dimensional spline with following properties :
 * KnotVector1 = {0, 0, 0, 1, 1, 1}
 * KnotVector2 = {0, 0, 0, 1, 1, 1}
 * ControlPoints = {{-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 1},
 *                  {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}}
*/
class A2DBSpline : public Test {
 public:
  A2DBSpline() :
      parameter_space(std::make_shared<NiceMock<Mock2dParameterSpace>>()),
      physical_space(std::make_shared<NiceMock<Mock2dPhysicalSpace>>()) {
    spl::BSplineGenerator<2> b_spline_generator(physical_space, parameter_space);
    b_spline = std::make_unique<spl::BSpline<2>>(b_spline_generator);
  }

 protected:
  std::unique_ptr<spl::BSpline<2>> b_spline;
  std::shared_ptr<NiceMock<Mock2dParameterSpace>> parameter_space;
  std::shared_ptr<NiceMock<Mock2dPhysicalSpace>> physical_space;
};

TEST_F(A2DBSpline, Corner) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {0})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {1})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.0}, ParamCoord{0.33333}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim1) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {0})[0], -0.33333, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {1})[0], -1.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.33333}, ParamCoord{0.0}}, {2})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, Center) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.5}, ParamCoord{0.5}}, {2})[0], 0.25, 0.00005);
}

TEST_F(A2DBSpline, Random) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {0})[0], 0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {1})[0], -0.5, 0.00005);
  ASSERT_NEAR(b_spline->Evaluate({ParamCoord{0.75}, ParamCoord{0.25}}, {2})[0], 0.14063, 0.00005);
}

TEST_F(A2DBSpline, CornerDer10) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {2}, {1, 0})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, CornerDer01) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {0}, {0, 1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {1}, {0, 1})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.0}}, {2}, {0, 1})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der10) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {2}, {1, 0})[0], 0.888889, 0.00005);
}

TEST_F(A2DBSpline, EdgeDim0Der01) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {0}, {0, 1})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {1}, {0, 1})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.0}, ParamCoord{0.33333}}, {2}, {0, 1})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, CenterDer10) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {0}, {1, 0})[0], 2.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {1}, {1, 0})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.5}, ParamCoord{0.5}}, {2}, {1, 0})[0], 0.0, 0.00005);
}

TEST_F(A2DBSpline, RandomDer10) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {1, 0})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {1, 0})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {1, 0})[0], -0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer01) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {0, 1})[0], 0.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {0, 1})[0], 2.000, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {0, 1})[0], 0.375, 0.00005);
}

TEST_F(A2DBSpline, RandomDer12) { // NOLINT
  mock_2dphysicalSpace(physical_space);
  mock_2dparameterSpace(parameter_space);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {0}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {1}, {1, 2})[0], 0.0, 0.00005);
  ASSERT_NEAR(b_spline->EvaluateDerivative({ParamCoord{0.75}, ParamCoord{0.25}}, {2}, {1, 2})[0], 4.0, 0.00005);
}
