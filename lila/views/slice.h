#pragma once

#include <lila/common.h>
#include <utility>

namespace lila {

struct Slice {
  size_type begin;
  size_type end;
  size_type step = 1;
};

inline constexpr size_type END = -1;
inline constexpr Slice ALL = {0, END, 1};

inline std::pair<Slice, size_type> adjusted_slice_length(Slice slice,
                                                         size_type length) {

  assert(slice.step != 0);

  // Finding out the correct begin
  if (slice.begin < 0) {
    slice.begin += length;
    if (slice.begin < 0) {
      slice.begin = (slice.step < 0) ? END : 0;
    }
  } else if (slice.begin >= length) {
    slice.begin = (slice.step < 0) ? length - 1 : length;
  }

  // Finding out the correct begin
  if (slice.end < 0) {
    slice.end += length;
    if (slice.end < 0) {
      slice.end = (slice.step < 0) ? END : 0;
    }
  } else if (slice.end >= length) {
    slice.end = (slice.step < 0) ? length - 1 : length;
  }

  // finding out the new length
  if (slice.step < 0) {
    if (slice.end < slice.begin) {
      return {slice, (slice.begin - slice.end - 1) / (-slice.step) + 1};
    }
  } else {
    if (slice.begin < slice.end) {
      return {slice, (slice.end - slice.begin - 1) / slice.step + 1};
    }
  }
  return {slice, 0};
}

} // namespace lila
