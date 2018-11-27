/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "nurbs_generator.h"
#include "surface_generator.h"
#include "nurbs.h"

#include "gmock/gmock.h"

#include "numeric_settings.h"

using testing::Test;
using testing::DoubleEq;
using testing::DoubleNear;

class ASurface : public Test {
 public:
  ASurface() :
        degree1_{Degree{1}},
        degree2_{Degree{2}},
        knot_vector1_{
          std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0},
                                                             ParamCoord{1}, ParamCoord{1}}))},
        knot_vector2_{
          std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                             ParamCoord{0.25}, ParamCoord{0.25},
                                                             ParamCoord{0.5}, ParamCoord{0.5},
                                                             ParamCoord{0.75}, ParamCoord{0.75},
                                                             ParamCoord{1}, ParamCoord{1}, ParamCoord{1}}))},
        parameter_space1(std::make_shared<spl::ParameterSpace<1>>(spl::ParameterSpace<1>(knot_vector1_, degree1_))),
        parameter_space2(std::make_shared<spl::ParameterSpace<1>>(spl::ParameterSpace<1>(knot_vector2_, degree2_))) {
      control_points1 = {
          baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.0})),
          baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0}))
      };
      control_points2 = {
          baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
          baf::ControlPoint(std::vector<double>({0.0, 1.0, 1.0})),
          baf::ControlPoint(std::vector<double>({0.0, 0.0, 1.0})),
          baf::ControlPoint(std::vector<double>({0.0, -1.0, 1.0})),
          baf::ControlPoint(std::vector<double>({0.0, -1.0, 0.0})),
          baf::ControlPoint(std::vector<double>({0.0, -1.0, -1.0})),
          baf::ControlPoint(std::vector<double>({0.0, 0.0, -1.0})),
          baf::ControlPoint(std::vector<double>({0.0, 1.0, -1.0})),
          baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
      };
      weights1_ = {1.0, 1.0};
      weights2_ = {1, 0.70710, 1, 0.70710, 1, 0.70710, 1, 0.70710, 1};
      w_physical_space1 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points1,
                                                                                                        weights1_,
                                                                                                        {2}));
      w_physical_space2 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points2,
                                                                                                        weights2_,
                                                                                                        {9}));

      spl::NURBSGenerator<1> nurbs_generator1(w_physical_space1, parameter_space1);
      spl::NURBSGenerator<1> nurbs_generator2(w_physical_space2, parameter_space2);
      nurbs1 = std::make_shared<spl::NURBS<1>>(nurbs_generator1);
      nurbs2 = std::make_shared<spl::NURBS<1>>(nurbs_generator2);

      spl::SurfaceGenerator surfaceGenerator = spl::SurfaceGenerator(nurbs1, nurbs2);
      nurbsJoined = std::make_unique<spl::NURBS<2>>(surfaceGenerator);
  }

 protected:
  std::array<Degree, 1> degree1_;
  std::array<Degree, 1> degree2_;
  KnotVectors<1> knot_vector1_;
  KnotVectors<1> knot_vector2_;
  std::shared_ptr<spl::ParameterSpace<1>> parameter_space1;
  std::shared_ptr<spl::ParameterSpace<1>> parameter_space2;
  std::vector<baf::ControlPoint> control_points1;
  std::vector<baf::ControlPoint> control_points2;
  std::vector<double> weights1_;
  std::vector<double> weights2_;
  std::shared_ptr<spl::WeightedPhysicalSpace<1>> w_physical_space1;
  std::shared_ptr<spl::WeightedPhysicalSpace<1>> w_physical_space2;
  std::shared_ptr<spl::NURBS<1>> nurbs1;
  std::shared_ptr<spl::NURBS<1>> nurbs2;
  std::unique_ptr<spl::NURBS<2>> nurbsJoined;
};

TEST_F(ASurface, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(nurbs1->GetDimension(), 3);
  ASSERT_THAT(nurbs2->GetDimension(), 3);
  ASSERT_THAT(nurbsJoined->GetDimension(), 3);
}

TEST_F(ASurface, ReturnsCorrectKnotVector) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetKnotVector(1)->GetNumberOfKnots(), 12);
  ASSERT_THAT(nurbsJoined->GetKnotVector(0)->GetNumberOfKnots(), 4);
}

TEST_F(ASurface, ReturnCorrectNumberOfControlPoints) { //NOLINT
  ASSERT_THAT(nurbsJoined->GetNumberOfControlPoints(), 18);
}

TEST_F(ASurface, ReturnsCorrectWeights) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetWeight(std::array<int, 2>({1, 4})), DoubleEq(1.0));
  ASSERT_THAT(nurbsJoined->GetWeight(std::array<int, 2>({0, 3})), DoubleNear(0.70710, 0.0000001));
}

TEST_F(ASurface, ReturnCorrectControlPoints) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({1, 4}), 1), DoubleEq(-1));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({0, 3}), 2), DoubleEq(1));
}
