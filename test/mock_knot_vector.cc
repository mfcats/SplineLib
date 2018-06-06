/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "gmock/gmock.h"

#include "knot_vector.h"
#include "b_spline.h"

using testing::Test;
using testing::Eq;
using testing::DoubleEq;
using testing::Return;

class MockKnotVector : public baf::KnotVector {

  /* This mock knot vector is used to replace the knot vector {0, 0, 0, 1, 1, 1}.
   * The knot vector has the following methods:
   *    - int GetKnotSpan(ParamCoord param_coord)
   *    - bool IsLastKnot(ParamCoord param_coord)
   *    - bool IsInKnotVectorRange(ParamCoord param_coord)
   *    - ParamCoord knot(int index)
   *    - int NumberOfKnots()
   *    - ParamCoord GetLastKnot()
   *    - ConstKnotIterator begin()
   *    - ConstKnotIterator end()
   * */

 public:
  MOCK_CONST_METHOD1(GetKnotSpan, int(ParamCoord
      param_coord));
  MOCK_CONST_METHOD1(IsLastKnot, bool(ParamCoord
      param_coord));
  MOCK_CONST_METHOD1(IsInKnotVectorRange, bool(ParamCoord
      param_coord));
  MOCK_CONST_METHOD1(knot, ParamCoord(int
      index));
  MOCK_CONST_METHOD0(NumberOfKnots, int());
  MOCK_CONST_METHOD0(GetLastKnot, ParamCoord());
  MOCK_CONST_METHOD0(begin, baf::ConstKnotIterator());
  MOCK_CONST_METHOD0(end, baf::ConstKnotIterator());
};

class ABSplineWithMockKnotVector : public Test {
 public:
  ABSplineWithMockKnotVector() {
    std::array<baf::KnotVector, 1> knot_vector =
        {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                          ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
    std::array<int, 1> degree = {2};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.5, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.3})),
        baf::ControlPoint(std::vector<double>({3.0, 2.0})),
        baf::ControlPoint(std::vector<double>({4.0, 1.5})),
        baf::ControlPoint(std::vector<double>({4.0, 0.0}))
    };
    b_spline = std::make_unique<spl::BSpline<1>>(knot_vector, degree, control_points);
  }

 protected:
  std::unique_ptr<spl::BSpline<1>> b_spline;
};

TEST_F(ABSplineWithMockKnotVector, Returns0_0For0AndDim0) {
  MockKnotVector knotVector;
  EXPECT_CALL(knotVector, GetKnotSpan(ParamCoord{0.0}))
      .Times(1).WillOnce(Return(2));
  EXPECT_EQ(b_spline->Evaluate({ParamCoord{0.0}}, {0})[0], 0.0);
}