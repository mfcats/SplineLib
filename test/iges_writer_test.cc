/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include "gmock/gmock.h"
#include "iges_writer.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class AnIGESFileFromSpline : public Test {
 public:
  AnIGESFileFromSpline() : degree_{Degree{2}} {
      std::array<baf::KnotVector, 1> knot_vector =
          {baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{2}, ParamCoord{3},
                            ParamCoord{4}, ParamCoord{4}, ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})};
      control_points_ = {
          baf::ControlPoint(std::vector<double>({0.5, 1.0})),
          baf::ControlPoint(std::vector<double>({0.0, 1.0})),
          baf::ControlPoint(std::vector<double>({2.0, 9.0})),
          baf::ControlPoint(std::vector<double>({1.5, 7.5})),
          baf::ControlPoint(std::vector<double>({4.0, 2.3})),
          baf::ControlPoint(std::vector<double>({3.0, 3.0})),
          baf::ControlPoint(std::vector<double>({6.0, 2.5})),
          baf::ControlPoint(std::vector<double>({5.0, 0.0}))
      };
      knot_vector_[0] = {std::make_shared<baf::KnotVector>(knot_vector[0])};
      spl::BSplineGenerator<1> b_spline_generator(knot_vector_, degree_, control_points_);
      b_spline = std::make_unique<spl::BSpline<1>>(b_spline_generator);
      iges_writer_ = std::make_unique<io::IGESWriter>();
    }

    protected:
    std::unique_ptr<spl::BSpline<1>> b_spline;
    std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector_;
    std::array<Degree, 1> degree_;
    std::vector<baf::ControlPoint> control_points_;
    std::unique_ptr<io::IGESWriter> iges_writer_;
  };

TEST_F(AnIGESFileFromSpline, Test1) {
  iges_writer_->WriteIGESFile(*(b_spline), iges_write);
}
