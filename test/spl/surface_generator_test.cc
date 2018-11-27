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
#include "vtk_writer.h"

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
          baf::ControlPoint(std::vector<double>({10.0, 0.0, 0.0}))
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
      nbInter = 11;
      nbInterCmp = 2;
      spl::NURBSGenerator<1> nurbs_generator1(w_physical_space1, parameter_space1);
      spl::NURBSGenerator<1> nurbs_generator2(w_physical_space2, parameter_space2);
      nurbs1 = std::make_shared<spl::NURBS<1>>(nurbs_generator1);
      nurbs2 = std::make_shared<spl::NURBS<1>>(nurbs_generator2);

      spl::SurfaceGenerator surfaceGenerator = spl::SurfaceGenerator(nurbs1, nurbs2);
      nurbsJoined = std::make_shared<spl::NURBS<2>>(surfaceGenerator);
      std::vector<std::array<double, 3>> scaling = {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0},
          {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5},
          {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
      std::vector<std::array<double, 3>> scalingCmp(nbInterCmp, {1.0, 1.0, 1.0});

    spl::SurfaceGenerator surfaceGeneratorScaled = spl::SurfaceGenerator(nurbs1, nurbs2, nbInter, scaling);
    spl::SurfaceGenerator surfaceGeneratorScaledCmp = spl::SurfaceGenerator(nurbs1, nurbs2, nbInterCmp, scalingCmp);
    nurbsJoinedScaled = std::make_shared<spl::NURBS<2>>(surfaceGeneratorScaled);
    nurbsJoinedScaledCmp = std::make_shared<spl::NURBS<2>>(surfaceGeneratorScaledCmp);
  }

 protected:
  int nbInter, nbInterCmp;
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
  std::shared_ptr<spl::NURBS<2>> nurbsJoined;
  std::shared_ptr<spl::NURBS<2>> nurbsJoinedScaled;
  std::shared_ptr<spl::NURBS<2>> nurbsJoinedScaledCmp;
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

TEST_F(ASurface, ReturnsCorrectDiemnsionScaled) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetDimension(), 3);
}

TEST_F(ASurface, ReturnsCorrectNumberOfKnotsFirstDimension) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetKnotVector(0)->GetNumberOfKnots(), 13);
}

TEST_F(ASurface, ReturnsCorrectNumberOfControlPoints) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetNumberOfControlPoints(), 99);
}

TEST_F(ASurface, RetursCorrectControlPoint) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({1, 4}), 1), DoubleEq(-1));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({5, 4}), 1), DoubleEq(-0.5));
}

TEST_F(ASurface, ReturnCorrectControlPoint_0Dim) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({1, 4}), 0), DoubleEq(1));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({5, 4}), 0), DoubleEq(5));
}

TEST_F(ASurface, CompareBothMethods_nbControlPoints) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetNumberOfControlPoints(), nurbsJoinedScaledCmp->GetNumberOfControlPoints());
}

TEST_F(ASurface, CompareBothMethods_nbKnots) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetKnotVector(0)->GetNumberOfKnots(),
      nurbsJoinedScaledCmp->GetKnotVector(0)->GetNumberOfKnots());
}

TEST_F(ASurface, CompareBothMethods_ControlPoint) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({1, 6}), 0), DoubleEq(
      nurbsJoinedScaledCmp->GetControlPoint(std::array<int, 2>({1, 6}), 0)));
}

TEST_F(ASurface, ReturnCorrectVTK ) { // NOLINT
  std::vector<std::vector<int>> scattering_({{30, 30}});
  std::unique_ptr<io::VTKWriter> vtk_writer_(std::make_unique<io::VTKWriter>());
  std::vector<std::any> splines;
  std::vector<std::any> splines2;
  std::any nurbsJoined_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbsJoined);
  std::any nurbsJoinedScaledCmp_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbsJoinedScaledCmp);
  splines = {nurbsJoined_any};
  splines2 = {nurbsJoinedScaledCmp_any};
  remove("splines.vtk");
  remove("splines.vkt");
  vtk_writer_->WriteFile(splines, "splines.vtk", scattering_);
  vtk_writer_->WriteFile(splines2, "splines2.vtk", scattering_);
}

class AComplexSurface : public Test {
 public:
  AComplexSurface() :
      degree1_{Degree{3}},
      degree2_{Degree{2}},
      knot_vector1_{
          std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                             ParamCoord{0.25}, ParamCoord{0.5},
                                                             ParamCoord{0.75}, ParamCoord{1}, ParamCoord{1},
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
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 2.0}))
    };
    control_points2 = {
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.0})),
    };
    weights1_ = {1, 0.5, 0.5, 0.5, 0.5, 0.5, 1};
    weights2_ = {1, 0.70710, 1, 0.70710, 1, 0.70710, 1, 0.70710, 1};
    w_physical_space1 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points1,
                                                                                                      weights1_,
                                                                                                      {7}));
    w_physical_space2 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points2,
                                                                                                      weights2_,
                                                                                                      {9}));

    spl::NURBSGenerator<1> nurbs_generator1(w_physical_space1, parameter_space1);
    spl::NURBSGenerator<1> nurbs_generator2(w_physical_space2, parameter_space2);
    nurbs1 = std::make_shared<spl::NURBS<1>>(nurbs_generator1);
    nurbs2 = std::make_shared<spl::NURBS<1>>(nurbs_generator2);
    nbInter = 101;
    std::vector<std::array<double, 3>> scaling(nbInter, {1.0, 1.0, 1.0});
    spl::SurfaceGenerator surfaceGenerator = spl::SurfaceGenerator(nurbs1, nurbs2, nbInter, scaling);
    nurbsJoined = std::make_shared<spl::NURBS<2>>(surfaceGenerator);
    double fct;
    for (int i = 50; i < 90; ++i) {
      fct = i < 70 ? (1 - 0.04 * (i - 50)) : (1 - 0.04 * (90 - i));
      scaling[i] = {fct, fct, fct};
    }
    spl::SurfaceGenerator surfaceGeneratorScaled = spl::SurfaceGenerator(nurbs1, nurbs2, nbInter, scaling);
    nurbsJoinedScaled = std::make_shared<spl::NURBS<2>>(surfaceGeneratorScaled);
  }

 protected:
  std::array<Degree, 1> degree1_;
  std::array<Degree, 1> degree2_;
  int nbInter;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector1_;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector2_;
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
  std::shared_ptr<spl::NURBS<2>> nurbsJoined;
  std::shared_ptr<spl::NURBS<2>> nurbsJoinedScaled;
};

TEST_F(AComplexSurface, ReturnsCorrectKnotVectorSize) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetKnotVector(1)->GetNumberOfKnots(), 12);
  ASSERT_THAT(nurbsJoinedScaled->GetKnotVector(0)->GetNumberOfKnots(), nbInter + degree1_[0].get() + 1);
}

TEST_F(AComplexSurface, ReturnCorrectVTK ) { // NOLINT
  std::vector<std::vector<int>> scattering_({{100, 100}});
  std::unique_ptr<io::VTKWriter> vtk_writer_(std::make_unique<io::VTKWriter>());
  std::vector<std::any> splines;
  std::any nurbsJoinedScaled_any = std::make_any<std::shared_ptr<spl::NURBS<2>>>(nurbsJoinedScaled);
  splines = {nurbsJoinedScaled_any};
  remove("splines.vtk");
  vtk_writer_->WriteFile(splines, "splines.vtk", scattering_);
}

TEST_F(AComplexSurface , ReturnCorrectControlPoint_0_7) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({0, 0}), 0), DoubleEq(2.0));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({0, 0}), 1), DoubleEq(0.0));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({0, 0}), 2), DoubleEq(0.85));
}

TEST_F(AComplexSurface , ReturnCorrectControlPoint_4_5) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({4, 5}), 0), DoubleNear(1.71768, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({4, 5}), 1), DoubleNear(0.0915789, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({4, 5}), 2), DoubleNear(1.13107, 0.00001));
}

TEST_F(AComplexSurface, ReturnCorrectControlPoint_19_1) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({19, 1}), 0), DoubleNear(1.19567, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({19, 1}), 1), DoubleNear(0.0490178, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({19, 1}), 2), DoubleNear(0.4101, 0.00001));
}

TEST_F(AComplexSurface, ReturnCorrectControlPoint_77_1) { // NOLINT
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({77, 1}), 0), DoubleNear(0.277145, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({77, 1}), 1), DoubleNear(1.1, 0.00001));
  ASSERT_THAT(nurbsJoined->GetControlPoint(std::array<int, 2>({77, 1}), 2), DoubleNear(0.662169, 0.00001));
}

TEST_F(AComplexSurface, ReturnCorrectControlPointScaled_19_1) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({19, 1}), 0), DoubleNear(1.19567, 0.00001));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({19, 1}), 1), DoubleNear(0.0490178, 0.00001));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({19, 1}), 2), DoubleNear(0.4101, 0.00001));
}

TEST_F(AComplexSurface, ReturnCorrectControlPointScaled_77_1) { // NOLINT
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({77, 1}), 0), DoubleNear(0.201014, 0.00001));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({77, 1}), 1), DoubleNear(1.048, 0.00001));
  ASSERT_THAT(nurbsJoinedScaled->GetControlPoint(std::array<int, 2>({77, 1}), 2), DoubleNear(0.645199, 0.00001));
}

class AComplexSurface2 : public Test {
 public:
  AComplexSurface2() :
      degree1_{Degree{2}},
      degree2_{Degree{2}},
      knot_vector1_{
          std::make_shared<baf::KnotVector>(baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0},
                                                             ParamCoord{0.2}, ParamCoord{0.4},
                                                             ParamCoord{0.6}, ParamCoord{0.8}, ParamCoord{1},
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
        baf::ControlPoint(std::vector<double>({2.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 1.0})),
        baf::ControlPoint(std::vector<double>({1.0, 0.0, 0.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 1.0, 2.0})),
        baf::ControlPoint(std::vector<double>({1.0, 1.0, 2.0}))
    };
    control_points2 = {
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, 0.1})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, 0.0})),
        baf::ControlPoint(std::vector<double>({0.0, -0.15, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.0, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, -0.1})),
        baf::ControlPoint(std::vector<double>({0.0, 0.15, 0.0})),
    };
    weights1_ = {1, 0.5, 0.5, 0.5, 0.5, 0.5, 1};
    weights2_ = {1, 0.70710, 1, 0.70710, 1, 0.70710, 1, 0.70710, 1};
    w_physical_space1 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points1,
                                                                                                      weights1_,
                                                                                                      {7}));
    w_physical_space2 = std::make_shared<spl::WeightedPhysicalSpace<1>>(spl::WeightedPhysicalSpace<1>(control_points2,
                                                                                                      weights2_,
                                                                                                      {9}));

    spl::NURBSGenerator<1> nurbs_generator1(w_physical_space1, parameter_space1);
    spl::NURBSGenerator<1> nurbs_generator2(w_physical_space2, parameter_space2);
    nurbs1 = std::make_shared<spl::NURBS<1>>(nurbs_generator1);
    nurbs2 = std::make_shared<spl::NURBS<1>>(nurbs_generator2);
    nbInter = 101;
  }

 protected:
  std::array<Degree, 1> degree1_;
  std::array<Degree, 1> degree2_;
  int nbInter;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector1_;
  std::array<std::shared_ptr<baf::KnotVector>, 1> knot_vector2_;
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
};

TEST_F(AComplexSurface2, ThrowExceptionIfDegreeTooLow) { // NOLINT
  std::vector<std::array<double, 3>> scaling(nbInter, {1.0, 1.0, 1.0});
  ASSERT_THROW(spl::SurfaceGenerator surfaceGenerator =
      spl::SurfaceGenerator(nurbs1, nurbs2, nbInter, scaling), std::runtime_error);
}

