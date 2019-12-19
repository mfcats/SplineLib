/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <fstream>

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

TEST_F(A1DRandomBSplineToWriteForFailingTest, WritesCorrectXMLFile) {  // NOLINT
  ASSERT_THAT(testing::UnitTest::GetInstance()->failed_test_count(), Eq(0));

  io::XMLWriter xml_writer;
  xml_writer.WriteFile({std::make_any<std::shared_ptr<splinelib::src::spl::BSpline<1>>>(random_b_spline_)},
                       xml_file_with_random_spline);
  std::system(("cat " + std::string(xml_file_with_random_spline)).c_str());
  std::cout << std::endl;
  remove(xml_file_with_random_spline);

  testing::internal::CaptureStdout();
  random_spline_writer::WriteToXML<1>(random_b_spline_, testing::UnitTest::GetInstance());
  std::string output = testing::internal::GetCapturedStdout();
  std::string expected_message = "\nAt least one test in test case A1DRandomBSplineToWriteForFailingTest failed. Here "
                                 "is the tested B-spline in XML-format:\n\n";
  ASSERT_THAT(output.substr(0, expected_message.length()), expected_message);
  std::string xml_file_output = output.substr(expected_message.length(), output.length() - expected_message.length());
  std::cout << xml_file_output << std::endl;
  std::ofstream newFile;
  newFile.open("random.xml");
  newFile << xml_file_output;
  newFile.close();
  io::XMLReader xml_reader;
  auto written_b_spline = std::any_cast<std::shared_ptr<spl::BSpline<1>>>(xml_reader.ReadFile("random.xml")[0]);
  ASSERT_THAT(written_b_spline->AreEqual(*random_b_spline_), true);
  remove("random.xml");
}