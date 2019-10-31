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
#include "random_b_spline_generator.h"
#include "random_nurbs_generator.h"

using testing::Test;
using testing::DoubleNear;

using namespace splinelib::src;

class Random1DBSplineForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random1DBSplineForKnotInsertionAndRemoval()
      : removed_{0}, param_coord_(ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}) {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomBSplineGenerator<1> spline_generator(limits, 10, 3);
    spl::BSpline<1> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<1>>(b_spline);

    spl::BSpline<1> insertion_spline(b_spline);
    insertion_spline.InsertKnot(param_coord_, 0);
    after_insertion_ = std::make_shared<spl::BSpline<1>>(insertion_spline);

    spl::BSpline<1> removal_spline(insertion_spline);
    removed_ = removal_spline.RemoveKnot(param_coord_, 0, 1e-10);
    after_removal_ = std::make_shared<spl::BSpline<1>>(removal_spline);
  }

 protected:
  size_t removed_;
  ParametricCoordinate param_coord_;
  std::shared_ptr<spl::BSpline<1>> original_;
  std::shared_ptr<spl::BSpline<1>> after_insertion_;
  std::shared_ptr<spl::BSpline<1>> after_removal_;
};

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, HasOneMoreKnotAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 1);
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, HasOneMoreControlPointAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 1);
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 100; ++evaluated_point) {
    std::array<ParametricCoordinate, 1>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, HasRemovedKnot) {  // NOLINT
  ASSERT_THAT(removed_, 1);
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, HasOneLessKnotAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 1);
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, HasOneLessControlPointAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[0], after_insertion_->GetPointsPerDirection()[0] - 1);
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 100; ++evaluated_point) {
    std::array<ParametricCoordinate, 1>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random1DBSplineForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

class Random1DNURBSForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random1DNURBSForKnotInsertionAndRemoval()
      : removed_{0}, param_coord_(ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}) {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomNURBSGenerator<1> spline_generator(limits, 10, 3);
    spl::NURBS<1> nurbs(spline_generator);
    original_ = std::make_shared<spl::NURBS<1>>(nurbs);

    spl::NURBS<1> insertion_spline(nurbs);
    insertion_spline.InsertKnot(param_coord_, 0);
    after_insertion_ = std::make_shared<spl::NURBS<1>>(insertion_spline);

    spl::NURBS<1> removal_spline(insertion_spline);
    removed_ = removal_spline.RemoveKnot(param_coord_, 0, 1e-10);
    after_removal_ = std::make_shared<spl::NURBS<1>>(removal_spline);
  }

 protected:
  size_t removed_;
  ParametricCoordinate param_coord_;
  std::shared_ptr<spl::NURBS<1>> original_;
  std::shared_ptr<spl::NURBS<1>> after_insertion_;
  std::shared_ptr<spl::NURBS<1>> after_removal_;
};

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, HasOneMoreKnotAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 1);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, HasOneMoreControlPointAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 1);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 100; ++evaluated_point) {
    std::array<ParametricCoordinate, 1>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, HasRemovedKnot) {  // NOLINT
  ASSERT_THAT(removed_, 1);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, HasOneLessKnotAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 1);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, HasOneLessControlPointAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[0], after_insertion_->GetPointsPerDirection()[0] - 1);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

TEST_F(Random1DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 100; ++evaluated_point) {
    std::array<ParametricCoordinate, 1>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

class Random2DBSplineForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random2DBSplineForKnotInsertionAndRemoval() : removed_({0}), param_coord_{
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}} {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomBSplineGenerator<2> spline_generator(limits, 10, 3);
    spl::BSpline<2> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<2>>(b_spline);

    spl::BSpline<2> insertion_spline(b_spline);
    insertion_spline.InsertKnot(param_coord_[0], 0);
    insertion_spline.InsertKnot(param_coord_[1], 1);
    insertion_spline.InsertKnot(param_coord_[2], 1);
    after_insertion_ = std::make_shared<spl::BSpline<2>>(insertion_spline);

    spl::BSpline<2> removal_spline(insertion_spline);
    removed_[0] = removal_spline.RemoveKnot(param_coord_[0], 0, 1e-10);
    removed_[1] = removal_spline.RemoveKnot(param_coord_[1], 1, 1e-10);
    removed_[2] = removal_spline.RemoveKnot(param_coord_[2], 1, 1e-10);
    after_removal_ = std::make_shared<spl::BSpline<2>>(removal_spline);
  }

 protected:
  std::array<size_t, 3> removed_;
  std::array<ParametricCoordinate, 3> param_coord_;
  std::shared_ptr<spl::BSpline<2>> original_;
  std::shared_ptr<spl::BSpline<2>> after_insertion_;
  std::shared_ptr<spl::BSpline<2>> after_removal_;
};

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, HasOneMoreKnotPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(after_insertion_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots() + 2);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, HasOneMoreControlPointPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1] + 2);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 25; ++evaluated_point) {
    std::array<ParametricCoordinate, 2>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, HasRemovedKnots) {  // NOLINT
  ASSERT_THAT(removed_[0], 1);
  ASSERT_THAT(removed_[1], 1);
  ASSERT_THAT(removed_[2], 1);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, HasOneLessKnotPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(after_removal_->GetKnotVector(1)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(1)->GetNumberOfKnots() - 2);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, HasOneLessControlPointPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[0], after_insertion_->GetPointsPerDirection()[0] - 1);
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[1], after_insertion_->GetPointsPerDirection()[1] - 2);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

TEST_F(Random2DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 25; ++evaluated_point) {
    std::array<ParametricCoordinate, 2>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

class Random2DNURBSForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random2DNURBSForKnotInsertionAndRemoval() : removed_({0}), param_coord_{
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}} {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomNURBSGenerator<2> spline_generator(limits, 8, 3);
    spl::NURBS<2> nurbs(spline_generator);
    original_ = std::make_shared<spl::NURBS<2>>(nurbs);

    spl::NURBS<2> insertion_spline(nurbs);
    insertion_spline.InsertKnot(param_coord_[0], 0);
    insertion_spline.InsertKnot(param_coord_[1], 0);
    insertion_spline.InsertKnot(param_coord_[2], 1);
    after_insertion_ = std::make_shared<spl::NURBS<2>>(insertion_spline);

    spl::NURBS<2> removal_spline(insertion_spline);
    removed_[0] = removal_spline.RemoveKnot(param_coord_[0], 0, 1e-10);
    removed_[1] = removal_spline.RemoveKnot(param_coord_[1], 0, 1e-10);
    removed_[2] = removal_spline.RemoveKnot(param_coord_[2], 1, 1e-10);
    after_removal_ = std::make_shared<spl::NURBS<2>>(removal_spline);
  }

 protected:
  std::array<size_t, 3> removed_;
  std::array<ParametricCoordinate, 3> param_coord_;
  std::shared_ptr<spl::NURBS<2>> original_;
  std::shared_ptr<spl::NURBS<2>> after_insertion_;
  std::shared_ptr<spl::NURBS<2>> after_removal_;
};

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, HasOneMoreKnotPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 2);
  ASSERT_THAT(after_insertion_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots() + 1);
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, HasOneMoreControlPointPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 2);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1] + 1);
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 2>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, HasRemovedKnots) {  // NOLINT
  ASSERT_THAT(removed_[0], 1);
  ASSERT_THAT(removed_[1], 1);
  ASSERT_THAT(removed_[2], 1);
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, HasOneLessKnotPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 2);
  ASSERT_THAT(after_removal_->GetKnotVector(1)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(1)->GetNumberOfKnots() - 1);
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

TEST_F(Random2DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 2>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

class Random3DBSplineForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random3DBSplineForKnotInsertionAndRemoval() : removed_({0}), param_coord_{
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}} {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomBSplineGenerator<3> spline_generator(limits, 10, 3);
    spl::BSpline<3> b_spline(spline_generator);
    original_ = std::make_shared<spl::BSpline<3>>(b_spline);

    spl::BSpline<3> insertion_spline(b_spline);
    insertion_spline.InsertKnot(param_coord_[0], 0);
    insertion_spline.InsertKnot(param_coord_[1], 1);
    insertion_spline.InsertKnot(param_coord_[2], 2);
    after_insertion_ = std::make_shared<spl::BSpline<3>>(insertion_spline);

    spl::BSpline<3> removal_spline(insertion_spline);
    removed_[0] = removal_spline.RemoveKnot(param_coord_[0], 0, 1e-10);
    removed_[1] = removal_spline.RemoveKnot(param_coord_[1], 1, 1e-10);
    removed_[2] = removal_spline.RemoveKnot(param_coord_[2], 2, 1e-10);
    after_removal_ = std::make_shared<spl::BSpline<3>>(removal_spline);
  }

 protected:
  std::array<size_t, 3> removed_;
  std::array<ParametricCoordinate, 3> param_coord_;
  std::shared_ptr<spl::BSpline<3>> original_;
  std::shared_ptr<spl::BSpline<3>> after_insertion_;
  std::shared_ptr<spl::BSpline<3>> after_removal_;
};

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, HasOneMoreKnotPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(after_insertion_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(after_insertion_->GetKnotVector(2)->GetNumberOfKnots(),
              original_->GetKnotVector(2)->GetNumberOfKnots() + 1);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, HasOneMoreControlPointPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2] + 1);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 3>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, HasRemovedKnots) {  // NOLINT
  ASSERT_THAT(removed_[0], 1);
  ASSERT_THAT(removed_[1], 1);
  ASSERT_THAT(removed_[2], 1);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, HasOneLessKnotPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(after_removal_->GetKnotVector(1)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(1)->GetNumberOfKnots() - 1);
  ASSERT_THAT(after_removal_->GetKnotVector(2)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(2)->GetNumberOfKnots() - 1);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, HasOneLessControlPointPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[0], after_insertion_->GetPointsPerDirection()[0] - 1);
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[1], after_insertion_->GetPointsPerDirection()[1] - 1);
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[2], after_insertion_->GetPointsPerDirection()[2] - 1);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

TEST_F(Random3DBSplineForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 3>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

class Random3DNURBSForKnotInsertionAndRemoval : public Test {  // NOLINT
 public:
  Random3DNURBSForKnotInsertionAndRemoval() : removed_({0}), param_coord_{
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)},
      ParametricCoordinate{util::random::GetUniformRandom<double>(0.01, 0.99)}} {
    std::array<ParametricCoordinate, 2> limits = {ParametricCoordinate{0.0}, ParametricCoordinate{1.0}};
    spl::RandomNURBSGenerator<3> spline_generator(limits, 5, 3);
    spl::NURBS<3> nurbs(spline_generator);
    original_ = std::make_shared<spl::NURBS<3>>(nurbs);

    spl::NURBS<3> insertion_spline(nurbs);
    insertion_spline.InsertKnot(param_coord_[0], 0);
    insertion_spline.InsertKnot(param_coord_[1], 1);
    insertion_spline.InsertKnot(param_coord_[2], 2);
    after_insertion_ = std::make_shared<spl::NURBS<3>>(insertion_spline);

    spl::NURBS<3> removal_spline(insertion_spline);
    removed_[0] = removal_spline.RemoveKnot(param_coord_[0], 0, 1e-10);
    removed_[1] = removal_spline.RemoveKnot(param_coord_[1], 1, 1e-10);
    removed_[2] = removal_spline.RemoveKnot(param_coord_[2], 2, 1e-10);
    after_removal_ = std::make_shared<spl::NURBS<3>>(removal_spline);
  }

 protected:
  std::array<size_t, 3> removed_;
  std::array<ParametricCoordinate, 3> param_coord_;
  std::shared_ptr<spl::NURBS<3>> original_;
  std::shared_ptr<spl::NURBS<3>> after_insertion_;
  std::shared_ptr<spl::NURBS<3>> after_removal_;
};

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, HasOneMoreKnotPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetKnotVector(0)->GetNumberOfKnots(),
              original_->GetKnotVector(0)->GetNumberOfKnots() + 1);
  ASSERT_THAT(after_insertion_->GetKnotVector(1)->GetNumberOfKnots(),
              original_->GetKnotVector(1)->GetNumberOfKnots() + 1);
  ASSERT_THAT(after_insertion_->GetKnotVector(2)->GetNumberOfKnots(),
              original_->GetKnotVector(2)->GetNumberOfKnots() + 1);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, HasOneMoreControlPointPerDirectionAfterKnotInsertion) {  // NOLINT
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[0], original_->GetPointsPerDirection()[0] + 1);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[1], original_->GetPointsPerDirection()[1] + 1);
  ASSERT_THAT(after_insertion_->GetPointsPerDirection()[2], original_->GetPointsPerDirection()[2] + 1);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotInsertion) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 3>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_insertion_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, HasRemovedKnots) {  // NOLINT
  ASSERT_THAT(removed_[0], 1);
  ASSERT_THAT(removed_[1], 1);
  ASSERT_THAT(removed_[2], 1);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, HasOneLessKnotPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetKnotVector(0)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(0)->GetNumberOfKnots() - 1);
  ASSERT_THAT(after_removal_->GetKnotVector(1)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(1)->GetNumberOfKnots() - 1);
  ASSERT_THAT(after_removal_->GetKnotVector(2)->GetNumberOfKnots(),
              after_insertion_->GetKnotVector(2)->GetNumberOfKnots() - 1);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, HasOneLessControlPointPerDirectionAfterKnotRemoval) {  // NOLINT
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[0], after_insertion_->GetPointsPerDirection()[0] - 1);
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[1], after_insertion_->GetPointsPerDirection()[1] - 1);
  ASSERT_THAT(after_removal_->GetPointsPerDirection()[2], after_insertion_->GetPointsPerDirection()[2] - 1);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, DoesNotChangeAfterKnotInsertionAndRemoval) {  // NOLINT
  ASSERT_THAT(original_->AreEqual((*after_removal_.get()), 1e-10), true);
}

TEST_F(Random3DNURBSForKnotInsertionAndRemoval, DoesNotChangeGeometricallyAfterKnotRemoval) {  // NOLINT
  for (int evaluated_point = 0; evaluated_point < 10; ++evaluated_point) {
    std::array<ParametricCoordinate, 3>
        param_coord{ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)},
                    ParametricCoordinate{util::random::GetUniformRandom<double>(0.0, 1.0)}};
    for (int j = 0; j < original_->GetPointDim(); ++j) {
      ASSERT_THAT(after_removal_->Evaluate(param_coord, {j})[0],
                  DoubleNear(original_->Evaluate(param_coord, {j})[0], 1e-10));
    }
  }
}
