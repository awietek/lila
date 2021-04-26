// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_LOGGER_H_
#define LILA_LOGGER_H_

#include <iostream>
#include <cstdlib>

#define FMT_HEADER_ONLY
#include <lila/external/fmt/format.h>

namespace lila {

  class Logger
  {
  public:
    Logger() : verbosity_(0) {};

    void set_verbosity(int verbosity) { verbosity_ = verbosity; }
    int verbosity() { return verbosity_; } 

    template <typename... Args>
    void out(const std::string &format, const Args &... args) {
      std::cout << fmt::format(format, args...);
    }

    template <typename... Args>
    void warn(const std::string &format, const Args &... args) {
      std::cout << fmt::format(format, args...) << std::flush;
    }

    template <typename... Args>
    void err(const std::string &format, const Args &... args) {
      std::cerr << fmt::format(format, args...) << std::flush;
      exit(EXIT_FAILURE);
    }

    template <typename... Args>
    void out(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cout << fmt::format(format, args...);
    }

    template <typename... Args>
    void warn(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cout << fmt::format(format, args...) << std::flush;
    }

    template <typename... Args>
    void err(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cerr << fmt::format(format, args...) << std::flush;
      exit(EXIT_FAILURE);
    }
    
  private:
    int verbosity_;

  };

  extern Logger logger;

}
#endif
