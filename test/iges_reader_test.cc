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
#include "iges_reader.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class AnIGESReader : public Test {
 public:
  AnIGESReader() {
  }
};

TEST_F(AnIGESReader, ReadNURBSCurveFromIGESFile) {
  util::IGESReader reader = util::IGESReader("/Users/christophsusen/SplineLib/test/test.iges");
  std::shared_ptr<spl::Spline<1>> spline = reader.ReadIGESFile(4);
  ASSERT_THAT(spline->Evaluate({ParamCoord{0.0}}, {0})[0], DoubleNear(-2.23308, 0.005));
}


