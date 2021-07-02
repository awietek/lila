#pragma once
#include <chrono>
#include <lila/utils/logger.h>
#include <string>

namespace lila {
using namespace std::chrono;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

auto inline rightnow() -> decltype(high_resolution_clock::now()) {
  return high_resolution_clock::now();
}

template <class Clock,
          class Duration = typename Clock::duration> // Uff, that's awful, Fuck C++
void timing(time_point<Clock, Duration> const &t0, 
	    time_point<Clock, Duration> const &t1, std::string msg = "",
            int verbosity = 1) {
  auto td = duration_cast<milliseconds>(t1 - t0).count();
  if (msg != "")
    Log.out(verbosity, "{}: {}.{} secs", msg, td / 1000, td % 1000);
  else
    Log.out(verbosity, "{}.{} secs", td / 1000, td % 1000);
}

inline void tic(bool begin = true, std::string msg = "", int verbosity = 1) {
  static auto t0 = rightnow();
  if (begin) {
    t0 = rightnow();
  } else {
    auto t1 = rightnow();
    timing(t0, t1, msg, verbosity);
  }
}

inline void toc(std::string msg = "", int verbosity = 1) {
  tic(false, msg, verbosity);
}

} // namespace lila
