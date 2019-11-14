/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetValue(CONTAINER_TYPE &container,  // NOLINT(runtime/references)
                                                         typename CONTAINER_TYPE::size_type index) {
#ifndef NDEBUG
  return container.at(index);
#else
  return container[index];
#endif
}

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const & GetValue(CONTAINER_TYPE const &container,
                                                               typename CONTAINER_TYPE::size_type index) {
#ifndef NDEBUG
  return container.at(index);
#else
  return container[index];
#endif
}

template<typename CONTAINER_TYPE, typename NAME>
constexpr typename CONTAINER_TYPE::value_type &GetValue(CONTAINER_TYPE &container,  // NOLINT(runtime/references)
                                                        NamedType<int, NAME> index) {
#ifndef NDEBUG
  return container.at(index.Get());
#else
  return container[index.Get()];
#endif
}

template<typename CONTAINER_TYPE, typename NAME>
constexpr typename CONTAINER_TYPE::value_type const & GetValue(CONTAINER_TYPE const &container,
                                                               NamedType<int, NAME> index) {
#ifndef NDEBUG
  return container.at(index.Get());
#else
  return container[index.Get()];
#endif
}

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetFront(CONTAINER_TYPE &container) {  // NOLINT(runtime/references)
#ifndef NDEBUG
  return container.at(0);
#else
  return container.front();
#endif
}

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const &GetFront(CONTAINER_TYPE const &container) {
#ifndef NDEBUG
  return container.at(0);
#else
  return container.front();
#endif
}

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type &GetBack(CONTAINER_TYPE &container) {  // NOLINT(runtime/references)
#ifndef NDEBUG
  return container.at(container.size() - 1);
#else
  return container.back();
#endif
}

template<typename CONTAINER_TYPE>
constexpr typename CONTAINER_TYPE::value_type const &GetBack(CONTAINER_TYPE const &container) {
#ifndef NDEBUG
  return container.at(container.size() - 1);
#else
  return container.back();
#endif
}
