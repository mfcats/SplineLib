/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_UTIL_RANDOM_H_
#define SRC_UTIL_RANDOM_H_

#include <random>

namespace util {
class Random {
 public:
  template<class T>
  static T GetBinomialRandom(double min, double max, double distance) {
    static std::default_random_engine generator;
    std::binomial_distribution<int> distribution(static_cast<int>((max - min) / distance), 0.5);
    return distribution(generator) * distance + min;
  }

  template<class T>
  static T GetUniformRandom(double min, double max) {
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
  }
};
}  // namespace util

#endif  // SRC_UTIL_RANDOM_H_
