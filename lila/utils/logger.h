#pragma once

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
      std::cout << fmt::format(format, args...) << "\n" ;
    }

    template <typename... Args>
    void warn(const std::string &format, const Args &... args) {
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
    }

    template <typename... Args>
    void err(const std::string &format, const Args &... args) {
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
      exit(EXIT_FAILURE);
    }

    template <typename... Args>
    void out(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cout << fmt::format(format, args...) << "\n" ;
    }

    template <typename... Args>
    void warn(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cout << fmt::format(format, args...) << "\n" << std::flush;
    }

    template <typename... Args>
    void err(int level, const std::string &format, const Args &... args) {
      if (level <= verbosity_)
	std::cerr << fmt::format(format, args...) << "\n" << std::flush;
      exit(EXIT_FAILURE);
    }
    
  private:
    int verbosity_;

  };

  inline Logger Log;

}
