#pragma once

namespace lila {

  struct Slice {
    size_type begin, end;
  }

  struct SliceStep {
    size_type begin, end, step;
  }

}
