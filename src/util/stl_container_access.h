/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_UTIL_STL_CONTAINER_ACCESS_H_
#define SRC_UTIL_STL_CONTAINER_ACCESS_H_

#include "src/util/named_type.h"

namespace splinelib::src::util::stl_container_access {
template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetValue(CONTAINER_TYPE &container,  // NOLINT(runtime/references)
                                                        typename CONTAINER_TYPE::size_type index);
template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const & GetValue(CONTAINER_TYPE const &container,
                                                               typename CONTAINER_TYPE::size_type index);

template<typename CONTAINER_TYPE, typename NAME>
constexpr typename CONTAINER_TYPE::value_type &GetValue(CONTAINER_TYPE &container,  // NOLINT(runtime/references)
                                                        NamedType<int, NAME> index);
template<typename CONTAINER_TYPE, typename NAME>
constexpr typename CONTAINER_TYPE::value_type const & GetValue(CONTAINER_TYPE const &container,
                                                               NamedType<int, NAME> index);

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetFront(CONTAINER_TYPE &container);  // NOLINT(runtime/references)
template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const &GetFront(CONTAINER_TYPE const &container);
template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetBack(CONTAINER_TYPE &container);  // NOLINT(runtime/references)
template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const &GetBack(CONTAINER_TYPE const &container);

#include "src/util/stl_container_access.inc"
}  // namespace splinelib::src::util::stl_container_access

namespace splinelib::src {
using util::stl_container_access::GetValue;
using util::stl_container_access::GetFront;
using util::stl_container_access::GetBack;
}  // namespace splinelib::src

#endif  // SRC_UTIL_STL_CONTAINER_ACCESS_H_
