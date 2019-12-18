/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include "gmock/gmock.h"

#include "src/io/xml_reader.h"
#include "test/spl/random/random_spline_utils.h"
#include "test/spl/random/write_random_spline.h"

using testing::Test;
using testing::Eq;

using namespace splinelib::test;

class A1DRandomBSplineToWriteForFailingTest : public Test {  // NOLINT
 public:
  A1DRandomBSplineToWriteForFailingTest()
      : random_b_spline_(RandomSplineUtils<1>::GenerateRandomBSpline(limits_, max_degree_, dimension_)) {}

 protected:
  const std::array<ParametricCoordinate, 2> limits_{{ParametricCoordinate{0}, ParametricCoordinate{1}}};
  const int max_degree_{4};
  const int dimension_{2};
  std::shared_ptr<spl::BSpline<1>> random_b_spline_;
};

std::string GetCommandOutput(const std::string &command);

TEST_F(A1DRandomBSplineToWriteForFailingTest, WritesCorrectXMLFile) {  // NOLINT
  ASSERT_THAT(testing::UnitTest::GetInstance()->failed_test_count(), Eq(0));
  random_spline_writer::WriteToXML<1>(random_b_spline_, testing::UnitTest::GetInstance());
  std::unique_ptr<io::XMLReader> xml_reader = std::make_unique<io::XMLReader>();
  std::string directory = GetCommandOutput("pwd");
  directory.erase(directory.length() - 1);
  const char* path_to_xml_file =
      (directory + "/spl/A1DRandomBSplineToWriteForFailingTest_WritesCorrectXMLFile.xml").c_str();
  auto b_spline_after = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader->ReadFile(path_to_xml_file)[0]);
  ASSERT_THAT(b_spline_after->AreEqual(*random_b_spline_, 10e-10), true);
  remove(path_to_xml_file);
}