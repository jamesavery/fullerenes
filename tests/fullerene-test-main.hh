// fullerene-test-main.hh -- the entry point for test binaries that need a
// flag GoogleTest does not have.
//
//   --sizes N,N,...   the fullerene sizes a size-swept test runs over.  The
//                     compiled-in default must be a quick, complete sweep
//                     (CTest passes no arguments); this widens it for a
//                     manual run.  Tunables are CLI flags, never environment
//                     variables -- see the project rules.
//   --database PATH   the isomer-database root, i.e. IsomerDB::database_path.
//                     The configured default points into the source tree,
//                     which need not hold the data.
//   -h, --help        this, plus a pointer to the standard --gtest_* flags.
//
// Anything else is handed to GoogleTest, so --gtest_filter and friends keep
// working.  A binary's main is then one line:
//
//   int main(int argc, char** argv) { return fullerene_test::run(argc, argv, {20, 24, 26}); }
//
// This is the single definition: claude-projects/fullerene-claude-shared's
// subproject_test_main.hh includes this file, so the sub-project binaries and
// the library's own share one CLI.
#pragma once

#include <gtest/gtest.h>

#include "fullerenes/isomerdb.hh"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace fullerene_test {

// The fullerene sizes a size-parameterised test sweeps: set by run() from the
// binary's default or the --sizes override, read inside the tests.
inline std::vector<int>& sizes() { static std::vector<int> s; return s; }

inline std::vector<int> parse_sizes(const std::string& csv) {
  std::vector<int> out;
  for (std::size_t i = 0; i < csv.size();) {
    std::size_t j = csv.find(',', i);
    if (j == std::string::npos) j = csv.size();
    if (j > i) out.push_back(std::atoi(csv.substr(i, j - i).c_str()));
    i = j + 1;
  }
  return out;
}

inline void print_help(const char* prog, const char* project) {
  std::printf("%s -- %s test binary\n\n", prog, project);
  std::printf("usage: %s [--sizes N,N,...] [--database PATH] [--gtest_* ...]\n\n", prog);
  std::printf("options:\n");
  std::printf("  %-22s %s", "--sizes N,N,...", "fullerene sizes a size-swept test runs (default: ");
  for (std::size_t k = 0; k < sizes().size(); k++) std::printf("%s%d", k ? "," : "", sizes()[k]);
  std::printf(")\n");
  std::printf("  %-22s %s\n", "--database PATH", "isomer-database root (default: " );
  std::printf("  %-22s   %s)\n", "", IsomerDB::database_path.c_str());
  std::printf("  %-22s %s\n", "-h, --help", "show this help and exit");
  std::printf("  %-22s %s\n", "--gtest_*", "standard GoogleTest flags (e.g. --gtest_filter=Suite.Case)");
}

// Consume the flags above, hand the rest to GoogleTest, run.  An unrecognised
// --flag is rejected rather than passed through: GoogleTest ignores what it
// does not know, so a silently swallowed --database would run the sweep
// against a different corpus than the caller named and report it as a pass.
inline int run(int argc, char** argv, std::vector<int> default_sizes,
               const char* project = "fullerene") {
  sizes() = std::move(default_sizes);
  std::vector<char*> rest{argv[0]};
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help")        { print_help(argv[0], project); return 0; }
    if (a == "--sizes" && i + 1 < argc)    { sizes() = parse_sizes(argv[++i]); continue; }
    if (a.rfind("--sizes=", 0) == 0)       { sizes() = parse_sizes(a.substr(8)); continue; }
    if (a == "--database" && i + 1 < argc) { IsomerDB::database_path = argv[++i]; continue; }
    if (a.rfind("--database=", 0) == 0)    { IsomerDB::database_path = a.substr(11); continue; }
    if (a.rfind("--gtest", 0) != 0 && a.rfind("--", 0) == 0) {
      std::fprintf(stderr, "%s: unknown option %s\n", argv[0], a.c_str());
      print_help(argv[0], project);
      return 2;
    }
    rest.push_back(argv[i]);
  }
  rest.push_back(nullptr);   // GoogleTest reads argv[argc]
  int rc = int(rest.size()) - 1;
  ::testing::InitGoogleTest(&rc, rest.data());
  return RUN_ALL_TESTS();
}

}  // namespace fullerene_test
