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
#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "gmock/gmock.h"

#include "src/spl/spline.h"
#include "src/io/xml_writer.h"

// The functions below is used to write a random spline to an XML-file whenever a test fails.
namespace splinelib::test::random_spline_writer {
template<int PARAMETRIC_DIMENSIONALITY>
static void WriteToXML(std::shared_ptr<splinelib::src::spl::Spline<PARAMETRIC_DIMENSIONALITY>> spl,
                       testing::UnitTest* test_instance) {
  std::string filename = "spl/" + std::string(test_instance->current_test_case()->name()) + "_" +
                         std::string(test_instance->current_test_info()->name()) + ".xml";
  splinelib::src::io::XMLWriter xml_writer;
  xml_writer.WriteFile(
      {std::make_any<std::shared_ptr<splinelib::src::spl::Spline<PARAMETRIC_DIMENSIONALITY>>>(spl)}, filename.c_str());
  std::cout << std::endl << "Random spline that lead to failing test was written to following XML-file:" << std::endl;
  std::cout <<  "build/test/" + filename << std::endl;
}
}  // namespace splinelib::test::random_spline_writer

#endif  // TEST_SPL_RANDOM_WRITE_RANDOM_SPLINE_H_
