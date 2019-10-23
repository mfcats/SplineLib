/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

// NOLINT(build/header_guard)

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::MultiIndexHandler(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index_length) :
    overflowed_{false}, multi_index_length_(multi_index_length) {
  for (auto &integer : current_multi_index_value_) integer = 0;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::MultiIndexHandler(MultiIndexHandler &&other) noexcept :
    overflowed_{false}, multi_index_length_(std::move(other.multi_index_length_)),
    current_multi_index_value_(std::move(other.current_multi_index_value_)) {}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> & MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator=(
    MultiIndexHandler &&rhs) noexcept {
  overflowed_ = std::move(rhs.overflowed_);
  multi_index_length_ = std::move(rhs.multi_index_length_);
  current_multi_index_value_ = std::move(rhs.current_multi_index_value_);
  return (*this);
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::begin() {
  MultiIndexHandler initialized_multi_index_handler(multi_index_length_);
  return initialized_multi_index_handler;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::end() {
  MultiIndexHandler multi_index_handler_with_highest_indices(multi_index_length_);
  multi_index_handler_with_highest_indices.overflowed_ = true;
  return multi_index_handler_with_highest_indices;
}

template<int PARAMETRIC_DIMENSIONALITY>
inline int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator[](Dimension const &dimension) const {
  return GetValue(current_multi_index_value_, dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> & MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator++() {
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++current_dimension) {
    int &current_value = GetValue(current_multi_index_value_, current_dimension);
    int const &current_length = GetValue(multi_index_length_, current_dimension);
    if ((current_value == (current_length - 1)) || (current_length == 0)) {
      if (current_dimension == Dimension{PARAMETRIC_DIMENSIONALITY - 1}) overflowed_ = true;
      current_value = 0;
    } else {
      ++current_value;
      break;
    }
  }
  return (*this);
}

template<int PARAMETRIC_DIMENSIONALITY>
inline const MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator++(int) {
  MultiIndexHandler result(*this);
  ++(*this);
  return result;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> & MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator--() {
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++current_dimension) {
    int &current_value = GetValue(current_multi_index_value_, current_dimension);
    if (current_value == 0) {
      if (current_dimension == Dimension{PARAMETRIC_DIMENSIONALITY - 1}) overflowed_ = false;
      current_value = (GetValue(multi_index_length_, current_dimension) - 1);
    } else {
      --current_value;
      break;
    }
  }
  return (*this);
}

template<int PARAMETRIC_DIMENSIONALITY>
inline const MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator--(int) {
  MultiIndexHandler result(*this);
  --(*this);
  return result;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> & MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator+(int count) {
  for (int i = 0; i < count; ++i) ++(*this);
  return (*this);
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> & MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator-(int count) {
  for (int i = 0; i < count; ++i) --(*this);
  return (*this);
}

template<int PARAMETRIC_DIMENSIONALITY>
bool MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator==(const MultiIndexHandler &rhs) const {
  return this->overflowed_ == rhs.overflowed_ && this->multi_index_length_ == rhs.multi_index_length_ &&
      this->current_multi_index_value_ == rhs.current_multi_index_value_;
}

template<int PARAMETRIC_DIMENSIONALITY>
bool MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::operator!=(const MultiIndexHandler &rhs) const {
  return this->overflowed_ != rhs.overflowed_ || this->multi_index_length_ != rhs.multi_index_length_ ||
      this->current_multi_index_value_ != rhs.current_multi_index_value_;
}

template<int PARAMETRIC_DIMENSIONALITY>
inline std::array<int, PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetCurrentIndex()
    const {
  return current_multi_index_value_;
}

template<int PARAMETRIC_DIMENSIONALITY>
inline int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetCurrent1DIndex() const {
  return Get1DIndex(multi_index_length_, current_multi_index_value_);
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::Get1DIndex(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &length,
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) {
  int index_1d = 0;
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++current_dimension) {
    int const current_slice_size = GetCurrentSliceSize(length, current_dimension);
    index_1d += (GetValue(multi_index, current_dimension) * current_slice_size);
  }
  return index_1d;
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::Get1DIndex(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  return Get1DIndex(multi_index_length_, multi_index);
}

template<int PARAMETRIC_DIMENSIONALITY>
inline void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::SetCurrentIndex(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) {
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++current_dimension) {
    if (GetValue(multi_index, current_dimension) < 0)
      throw std::invalid_argument("All indices have to be non-negative.");
    if (GetValue(multi_index, current_dimension) >= GetValue(multi_index_length_, current_dimension))
      throw std::invalid_argument("Each index has to be smaller than the given corresponding maximum index.");
  }
  current_multi_index_value_ = multi_index;
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::SetCurrentIndex(int index_1d) {
  if (index_1d < 0) throw std::invalid_argument("The 1D index has to be non-negative.");
  if (index_1d >= GetNumberOfTotalMultiIndices())
    throw std::invalid_argument("The 1D index has to be smaller than the number of possible multi-indices.");
  for (Dimension current_dimension{PARAMETRIC_DIMENSIONALITY - 1}; current_dimension >= Dimension{0};
       --current_dimension) {
    int const current_slice_size = GetCurrentSliceSize(multi_index_length_, current_dimension);
    int &current_value = GetValue(current_multi_index_value_, current_dimension);
    current_value = static_cast<int>(index_1d / current_slice_size);
    index_1d -= (current_value * current_slice_size);
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetNumberOfTotalMultiIndices() const {
  Dimension const maximum_dimension{PARAMETRIC_DIMENSIONALITY - 1};
  return (GetValue(multi_index_length_, maximum_dimension) *
          GetCurrentSliceSize(multi_index_length_, maximum_dimension));
}

template<int PARAMETRIC_DIMENSIONALITY>
std::array<int, PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetComplementaryIndex() const {
  std::array<int, PARAMETRIC_DIMENSIONALITY> multi_index{};
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY}; ++current_dimension) {
    GetValue(multi_index, current_dimension) = (GetValue(multi_index_length_, current_dimension) -
                                                GetValue(current_multi_index_value_, current_dimension) - 1);
  }
  return multi_index;
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::CollapseDimension(Dimension const &dimension) const {
  std::array<int, PARAMETRIC_DIMENSIONALITY - 1> multi_index{}, length{};
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY - 1};
       ++current_dimension) {
    int &current_index = GetValue(multi_index, current_dimension);
    int &current_length = GetValue(length, current_dimension);
    if (current_dimension < dimension) {
      current_index = GetValue(current_multi_index_value_, current_dimension);
      current_length = GetValue(multi_index_length_, current_dimension);
    } else if (current_dimension >= dimension) {
      Dimension const next_dimension = (current_dimension + Dimension{1});
      current_index = GetValue(current_multi_index_value_, next_dimension);
      current_length = GetValue(multi_index_length_, next_dimension);
    }
  }
  return MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1>::Get1DIndex(length, multi_index);
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetCurrentSliceSize(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &length, Dimension const &maximum_dimension) {
  int current_slice_size = 1;
  for (Dimension dimension{0}; dimension < maximum_dimension; ++dimension) {
    current_slice_size *= GetValue(length, dimension);
  }
  return current_slice_size;
}
