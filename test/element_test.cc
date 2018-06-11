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

#include "element_generator.h"

using testing::Test;
using testing::DoubleEq;

class A1DElement : public Test {
 public:
  A1DElement() : element(1, {0.5, 1.0}) {}

 protected:
  elm::Element element;
};

TEST_F(A1DElement, ReturnsCorrectDimension) { // NOLINT
  ASSERT_THAT(element.dimension(), 1);
}

TEST_F(A1DElement, ReturnsCorrectNumberOfNodes) { // NOLINT
  ASSERT_THAT(element.numberOfNodes(), 2);
}

TEST_F(A1DElement, ReturnsCorrectNode) { // NOLINT
  ASSERT_THAT(element.node(0), DoubleEq(0.5));
  ASSERT_THAT(element.node(1), DoubleEq(1.0));
}

class A1DElementGenerator : public Test {
 public:
  A1DElementGenerator() : element_generator(2,
                                            baf::KnotVector({ParamCoord{0}, ParamCoord{0}, ParamCoord{0}, ParamCoord{1},
                                                             ParamCoord{2}, ParamCoord{3}, ParamCoord{4}, ParamCoord{4},
                                                             ParamCoord{5}, ParamCoord{5}, ParamCoord{5}})) {}

 protected:
  elm::ElementGenerator element_generator;
};

TEST_F(A1DElementGenerator, ReturnsCorrectNumberOfElements) { // NOLINT
  ASSERT_THAT(element_generator.GetElementList().size(), 5);
}

TEST_F(A1DElementGenerator, Returns1DElements) { // NOLINT
  for (auto &element : element_generator.GetElementList()) {
    ASSERT_THAT(element.dimension(), 1);
  }
}

TEST_F(A1DElementGenerator, ReturnsElementsWith2Nodes) { // NOLINT
  for (auto &element : element_generator.GetElementList()) {
    ASSERT_THAT(element.numberOfNodes(), 2);
  }
}

TEST_F(A1DElementGenerator, ReturnsElementsWithCorrectNodes) { // NOLINT
  auto element_list = element_generator.GetElementList();
  for (int element = 0; element < element_list.size(); element++) {
    ASSERT_THAT(element_list[element].node(0), element);
    ASSERT_THAT(element_list[element].node(1), element + 1);
  }
}
