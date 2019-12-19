/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef TEST_SPL_RANDOM_WRITE_RANDOM_SPLINE_H_
#define TEST_SPL_RANDOM_WRITE_RANDOM_SPLINE_H_

#include <any>
#include <config_random_xml.h>
#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "gmock/gmock.h"

#include "src/spl/spline.h"
#include "src/io/xml_writer.h"

// The functions below are used to write a random spline in XML-format to console whenever a test fails.
namespace splinelib::test::random_spline_writer {
template<int PARAMETRIC_DIMENSIONALITY>
static void WriteToXML(std::vector<std::shared_ptr<splinelib::src::spl::BSpline<PARAMETRIC_DIMENSIONALITY>>> splines,
                       testing::UnitTest* test_instance) {
  splinelib::src::io::XMLWriter xml_writer;
  std::vector<std::any> any_splines;
  any_splines.reserve(splines.size());
  for (int spline_index = 0; spline_index < static_cast<int>(splines.size()); ++spline_index) {
    any_splines.emplace_back(
        std::make_any<std::shared_ptr<splinelib::src::spl::BSpline<PARAMETRIC_DIMENSIONALITY>>>(splines[spline_index]));
  }
  xml_writer.WriteFile(any_splines, xml_file_with_random_spline);
  std::cout << std::endl << "At least one test in test case " << test_instance->current_test_case()->name() <<
               " failed. Here is the tested B-spline (and the modified version(s)) in XML-format:" << std::endl <<
               std::endl;
  std::system(("cat " + std::string(xml_file_with_random_spline)).c_str());
  std::cout << std::endl;
  remove(xml_file_with_random_spline);
}

template<int PARAMETRIC_DIMENSIONALITY>
static void WriteToXML(std::vector<std::shared_ptr<splinelib::src::spl::NURBS<PARAMETRIC_DIMENSIONALITY>>> splines,
                       testing::UnitTest* test_instance) {
  splinelib::src::io::XMLWriter xml_writer;
  std::vector<std::any> any_splines;
  any_splines.reserve(splines.size());
  for (int spline_index = 0; spline_index < static_cast<int>(splines.size()); ++spline_index) {
    any_splines.emplace_back(
        std::make_any<std::shared_ptr<splinelib::src::spl::NURBS<PARAMETRIC_DIMENSIONALITY>>>(splines[spline_index]));
  }
  xml_writer.WriteFile(any_splines, xml_file_with_random_spline);
  std::cout << std::endl << "At least one test in test case " << test_instance->current_test_case()->name() <<
               " failed. Here is the tested NURBS (and the modified version(s)) in XML-format:" << std::endl <<
               std::endl;
  std::system(("cat " + std::string(xml_file_with_random_spline)).c_str());
  std::cout << std::endl;
  remove(xml_file_with_random_spline);
}
}  // namespace splinelib::test::random_spline_writer

#endif  // TEST_SPL_RANDOM_WRITE_RANDOM_SPLINE_H_
