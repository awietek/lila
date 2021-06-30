#pragma once

#include <lila/common.h>

namespace lila {

struct Slice {
  size_type begin;
  size_type end;
  size_type step = 1;
};

} // namespace lila
