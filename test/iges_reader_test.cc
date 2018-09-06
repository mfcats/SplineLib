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
#include "spline_generator.h"
#include "b_spline_generator.h"
#include "iges_1d_bspline_generator.h"
#include "iges_2d_bspline_generator.h"
#include "iges_1d_nurbs_generator.h"
#include "iges_2d_nurbs_generator.h"
#include <config.h>

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class AnIGESReader : public Test {
 public:
  AnIGESReader() {
  }
};

TEST_F(AnIGESReader, Read1DBSplineFromIGESFile) { // NOLINT
  spl::IGES1DBSplineGenerator reader = spl::IGES1DBSplineGenerator(path_to_iges_file);
  reader.ReadIGESFile(4);
  std::unique_ptr<spl::BSpline<1>> spline = std::make_unique<spl::BSpline<1>>(reader);
  ASSERT_THAT(spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleNear(-2.23308, 0.0005));
  ASSERT_THAT(spline->Evaluate({ParamCoord{0.0}}, {1})[0], DoubleNear(-0.01433, 0.0005));
  ASSERT_THAT(spline->Evaluate({ParamCoord{0.0}}, {2})[0], DoubleNear(-0.51255, 0.0005));
  ASSERT_THAT(spline->Evaluate({ParamCoord{1.0}}, {0})[0], DoubleNear(-1.3353, 0.0005));
  ASSERT_THAT(spline->Evaluate({ParamCoord{1.0}}, {1})[0], DoubleNear(0.450443, 0.0005));
  ASSERT_THAT(spline->Evaluate({ParamCoord{1.0}}, {2})[0], DoubleNear(-0.023586, 0.0005));
}

TEST_F(AnIGESReader, Read2DNURBSFromIGESFile) { // NOLINT
  spl::IGES2DNURBSGenerator reader2 = spl::IGES2DNURBSGenerator(path_to_iges_file);
  reader2.ReadIGESFile(2);
  std::unique_ptr<spl::NURBS<2>> spline2 = std::make_unique<spl::NURBS<2>>(reader2);
  ASSERT_THAT(spline2->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {0})[0], DoubleNear(0, 0.0005));
  ASSERT_THAT(spline2->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {1})[0], DoubleNear(0, 0.0005));
  ASSERT_THAT(spline2->Evaluate({ParamCoord{0.0}, ParamCoord{0.0}}, {2})[0], DoubleNear(-1, 0.0005));
  ASSERT_THAT(spline2->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {0})[0], DoubleNear(0, 0.0005));
  ASSERT_THAT(spline2->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {1})[0], DoubleNear(0, 0.0005));
  ASSERT_THAT(spline2->Evaluate({ParamCoord{1.0}, ParamCoord{1.0}}, {2})[0], DoubleNear(0, 0.0005));
}


