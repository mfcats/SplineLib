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

#include "b_spline.h"
#include "nurbs.h"
#include "vtk_writer.h"

using testing::Test;

template<int DIM>
void PrintSpline(std::shared_ptr<spl::BSpline<DIM>> spline) {
  std::cout << std::endl << "degrees:" << std::endl;
  for (int i = 0; i < DIM; ++i) {
    std::cout << spline->GetDegree(i).get() << "   ";
  }
  std::cout << std::endl << std::endl << "knot vectors:" << std::endl;
  for (int i = 0; i < DIM; ++i) {
    auto kv = spline->GetKnotVector(i);
    for (const auto &knot : *kv) {
      std::cout << knot.get() << "  ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << "control points:" << std::endl;
  for (int i = 0; i < spline->GetNumberOfControlPoints(); ++i) {
    for (int j = 0; j < spline->GetPointDim(); ++j) {
      std::cout << spline->GetControlPoint({i}, j) << "  ";
    }
    std::cout << std::endl;
  }
}

class BSpline1DFig5_39 : public Test {  // NOLINT
 public:
  BSpline1DFig5_39() {
    std::array<Degree, 1> degree = {Degree{4}};
    KnotVectors<1> knot_vector_before = {std::make_shared<baf::KnotVector>(
        baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3},
                         ParamCoord{0.3}, ParamCoord{0.6}, ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1},
                         ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({-0.5, 0.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({-2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({-1.0, 1.5})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.5})),
        baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.5})),
        baf::ControlPoint(std::vector<double>({0.5, 0.0}))
    };
    bspline_1d_before_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
    bspline_1d_after_ = std::make_shared<spl::BSpline<1>>(knot_vector_before, degree, control_points);
  }

 protected:
  std::shared_ptr<spl::BSpline<1>> bspline_1d_before_;
  std::shared_ptr<spl::BSpline<1>> bspline_1d_after_;
};

TEST_F(BSpline1DFig5_39, ReducesDegreeFrom5To4Correctly) {  // NOLINT

  // Elevate the degree of the spline by one and then reduce it to obtain a spline that has not changed geometrically.
  bspline_1d_after_->ElevateDegree(0);
  bool successful = bspline_1d_after_->ReduceDegree(0, 0.0001);

  // Print degree, knot vector and control points of the degree reduced spline to the console.
  // PrintSpline(bspline_1d_after_);

  // Write spline after degree reduction to VTK file for visualization.
  /* std::vector<std::any> splines;
  splines.emplace_back(std::make_any<std::shared_ptr<spl::BSpline<1>>>(bspline_1d_after_));
  splines.emplace_back(std::make_any<std::shared_ptr<spl::BSpline<1>>>(bspline_1d_before_));
  io::VTKWriter vtk_writer;
  vtk_writer.WriteFile(splines, "/Users/christophsusen/Desktop/test.vtk", {{40}, {40}}); */

  // The geometry of the spline and its degree should have remained unchanged.
  ASSERT_THAT(successful, true);
  ASSERT_THAT(bspline_1d_after_->GetDegree(0).get(), bspline_1d_before_->GetDegree(0).get());
  ASSERT_THAT(bspline_1d_after_->AreGeometricallyEqual(*bspline_1d_before_), true);
}

class A2DBSplineForDegreeReduction : public Test {  // NOLINT
 public:
  A2DBSplineForDegreeReduction() {
    std::array<Degree, 2> degree = {Degree{2}, Degree{1}};
    KnotVectors<2> knot_vector_before = {
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0.3}, ParamCoord{0.6},
                             ParamCoord{0.6}, ParamCoord{1}, ParamCoord{1}, ParamCoord{1}})),
        std::make_shared<baf::KnotVector>(
            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{1}, ParamCoord{1}}))};
    std::vector<baf::ControlPoint> control_points = {
        baf::ControlPoint(std::vector<double>({0.0, 0.0})), baf::ControlPoint(std::vector<double>({0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({0.0, 2.0})), baf::ControlPoint(std::vector<double>({0.0, 3.0})),
        baf::ControlPoint(std::vector<double>({0.0, 4.0})), baf::ControlPoint(std::vector<double>({0.0, 5.0})),

        baf::ControlPoint(std::vector<double>({2.0, 0.0})), baf::ControlPoint(std::vector<double>({2.0, 1.0})),
        baf::ControlPoint(std::vector<double>({2.0, 2.0})), baf::ControlPoint(std::vector<double>({2.0, 3.0})),
        baf::ControlPoint(std::vector<double>({2.0, 4.0})), baf::ControlPoint(std::vector<double>({2.0, 5.0}))
    };
    original_ = std::make_shared<spl::BSpline<2>>(knot_vector_before, degree, control_points);
    after_reduction_ = std::make_shared<spl::BSpline<2>>(*original_);
  }

 protected:
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_reduction_;
};

TEST_F(A2DBSplineForDegreeReduction, ReducesDegreeFrom4To3Correctly) {  // NOLINT

  // Elevate the degree of the spline by one and then reduce it to obtain a spline that has not changed geometrically.
  // TODO: Degree reduction after just one degree elevation does not work.
  after_reduction_->ElevateDegree(0);
  after_reduction_->ElevateDegree(0);
  bool successful = after_reduction_->ReduceDegree(0, 0.001);

  // Write spline after degree reduction to VTK file for visualization.
  /* std::vector<std::any> splines;
  splines.emplace_back(std::make_any<std::shared_ptr<spl::BSpline<2>>>(original_));
  splines.emplace_back(std::make_any<std::shared_ptr<spl::BSpline<2>>>(after_reduction_));
  io::VTKWriter vtk_writer;
  vtk_writer.WriteFile(splines, "/Users/christophsusen/Desktop/test.vtk", {{40, 40}, {40, 40}}); */

  // The geometry of the spline and its degree should have remained unchanged.
  ASSERT_THAT(successful, true);
  ASSERT_THAT(after_reduction_->GetDegree(1).get(), original_->GetDegree(1).get());
  ASSERT_THAT(after_reduction_->AreGeometricallyEqual(*original_), true);
}
