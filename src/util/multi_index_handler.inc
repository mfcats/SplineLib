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
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index_length)
    : overflowed_{false}, multi_index_length_(multi_index_length) {
  for (auto &integer : current_multi_index_value_) integer = 0;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::MultiIndexHandler(MultiIndexHandler &&other) noexcept
    : overflowed_{false}, multi_index_length_(std::move(other.multi_index_length_)),
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
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::SubtractWithConstantDimension(Dimension const &dimension) {
  ThrowIfDimensionInvalid();
  MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1> multi_index_handler_without_constant_dimension =
      GetMultiIndexHandlerWithCollapsedDimension(dimension);
  BackslideCollapsedDimension(--multi_index_handler_without_constant_dimension, dimension);
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
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::begin() {
  MultiIndexHandler initialized_multi_index_handler(multi_index_length_);
  return initialized_multi_index_handler;
}

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::end() {
  MultiIndexHandler overflowed_multi_index_handler(multi_index_length_);
  overflowed_multi_index_handler.overflowed_ = true;
  return overflowed_multi_index_handler;
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
    ThrowIfIndexInvalid(GetValue(multi_index, current_dimension), current_dimension);
  }
  current_multi_index_value_ = multi_index;
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::SetCurrentIndex(int index_1d) {
  ThrowIf1DIndexInvalid(index_1d);
  for (Dimension current_dimension{PARAMETRIC_DIMENSIONALITY - 1}; current_dimension >= Dimension{0};
       --current_dimension) {
    int const current_slice_size = GetCurrentSliceSize(multi_index_length_, current_dimension);
    int &current_value = GetValue(current_multi_index_value_, current_dimension);
    current_value = static_cast<int>(index_1d / current_slice_size);
    index_1d -= (current_value * current_slice_size);
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::SetCurrentIndexForDimension(int value, Dimension const &dimension) {
  ThrowIfDimensionInvalid();
  ThrowIfIndexInvalid();
  GetValue(current_multi_index_value_, dimension) = value;
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
  MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1> multi_index_handler_with_collapsed_dimension =
      GetMultiIndexHandlerWithCollapsedDimension(dimension);
  return multi_index_handler_with_collapsed_dimension.GetCurrent1DIndex();
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetNumberOfMultiIndicesForCollapsedDimension(
    Dimension const &dimension) const {
  if (PARAMETRIC_DIMENSIONALITY == 0 || (PARAMETRIC_DIMENSIONALITY == 1 && dimension.Get() == 0))
    return PARAMETRIC_DIMENSIONALITY;
  MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1> multi_index_handler_with_collapsed_dimension =
      GetMultiIndexHandlerWithCollapsedDimension(dimension);
  return multi_index_handler_with_collapsed_dimension.GetNumberOfTotalMultiIndices();
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetCurrentSliceSize(Dimension const &maximum_dimension) const {
  return GetCurrentSliceSize(multi_index_length_, maximum_dimension);
}

template<int PARAMETRIC_DIMENSIONALITY>
int MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetCurrentSliceComplementFill(
    Dimension const &maximum_dimension) const {
  int current_slice_fill = 1;
  for (Dimension current_dimension{0}; current_dimension < maximum_dimension; ++current_dimension) {
    int const current_slice_size = GetCurrentSliceSize(current_dimension);
    current_slice_fill += (GetValue(GetComplementaryIndex(), current_dimension) * current_slice_size);
  }
  return current_slice_fill;
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

template<int PARAMETRIC_DIMENSIONALITY>
MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1>
    MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::GetMultiIndexHandlerWithCollapsedDimension(
        Dimension const &dimension) const {
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
  MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1> multi_index_handler_with_collapsed_dimension(length);
  multi_index_handler_with_collapsed_dimension.SetCurrentIndex(multi_index);
  return multi_index_handler_with_collapsed_dimension;
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::BackslideCollapsedDimension(
    MultiIndexHandler<PARAMETRIC_DIMENSIONALITY - 1> const &multi_index_handler_without_constant_dimension,
    Dimension const &dimension) {
  std::array<int, PARAMETRIC_DIMENSIONALITY - 1> multi_index =
      multi_index_handler_without_constant_dimension.GetCurrentIndex();
  for (Dimension current_dimension{0}; current_dimension < Dimension{PARAMETRIC_DIMENSIONALITY};
       ++current_dimension) {
    int &current_index = GetValue(current_multi_index_value_, current_dimension);
    if (current_dimension < dimension) {
      current_index = GetValue(multi_index, current_dimension);
    } else if (current_dimension > dimension) {
      Dimension const lower_dimension = (current_dimension - Dimension{1});
      current_index = GetValue(multi_index, lower_dimension);
    }
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::ThrowIfDimensionInvalid(Dimension const &dimension) const {
  if (dimension < Dimension{0} || dimension >= Dimension{PARAMETRIC_DIMENSIONALITY})
    throw std::invalid_argument("splinelib::src::util::MultiIndexHandler::ThrowIfDimensionInvalid: The given dimension "
                                "has to be between 0 and (PARAMETRIC_DIMENSIONALITY - 1).");
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::ThrowIfIndexInvalid(int index, Dimension const &dimension) const {
  if (index < 0 || index >= GetValue(multi_index_length_, dimension))
    throw std::invalid_argument("splinelib::src::util::MultiIndexHandler::ThrowIfIndexInvalid: The index has to "
                                "be non-negative and smaller than the multi-index length in the given dimension.");
}

template<int PARAMETRIC_DIMENSIONALITY>
void MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>::ThrowIf1DIndexInvalid(int index_1d) const {
  if (index_1d < 0 || index_1d >= GetNumberOfTotalMultiIndices())
    throw std::invalid_argument("splinelib::src::util::MultiIndexHandler::ThrowIf1DIndexInvalid: The 1D index has "
                                "to be non-negative and smaller than the number of possible multi-indices.");
}
