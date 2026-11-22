// stratonix command-line program. A thin layer over the library: every command calls the same
// functions the C API calls.
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cli_commands.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/package/pipeline.hpp"

using namespace stratonix;

namespace {

struct Command {
  std::string summary;
  std::string usage;
  std::function<int(cli::Args&)> run;
};

std::map<std::string, Command>& commands() {
  static std::map<std::string, Command> c = {
      {"validate", {"Check a model description against the schema and the semantic rules",
                    "stratonix validate <description.yaml>", cli::cmd_validate}},
      {"build", {"Build, check and export a model; write figures and a run manifest",
                 "stratonix build <description.yaml> -o <new folder> [--threads N] [--reference-kernels] "
                 "[--no-figures] [--formats a,b,...]", cli::cmd_build}},
      {"section", {"Draw a section along a polyline through a built model",
                   "stratonix section <description.yaml> -o <new folder> (--name <section> | --line x1,y1;x2,y2;...)",
                   cli::cmd_section}},
      {"state", {"Build the model at an event of its timeline (stage 2E) and write its package",
                 "stratonix state <description.yaml> --age <Ma> -o <new folder> [--no-figures] [--formats a,b,...]", cli::cmd_state}},
      {"export", {"Export a model in chosen formats",
                  "stratonix export <description.yaml> -o <new folder> [--formats petrel|irap_ascii,gocad_tsurf,...] [--grid-spacing m]",
                  cli::cmd_export}},
      {"package", {"Write a self-contained, checksummed delivery package",
                   "stratonix package <description.yaml> -o <new folder> [--config package.yaml] [--zip] "
                   "[--formats petrel,vtk,...] [--no-figures]", cli::cmd_package}},
      {"ensemble", {"Run the seeded ensemble of a description's uncertainty block (stage 3C)",
                    "stratonix ensemble <description.yaml> -o <new folder> [--members N] [--seed S] [--threads N] [--no-figures]", cli::cmd_ensemble}},
      {"fold-frame", {"Fit a fold frame to bedding orientations and predict bedding normals (stage 4A)",
                      "stratonix fold-frame <bedding.csv> [--foliation <axial_foliation.csv>] [--at <points.csv>] -o <frame.json>",
                      cli::cmd_fold_frame}},
      {"plugins", {"List the registered plugins", "stratonix plugins", cli::cmd_plugins}},
      {"bench", {"Run timing benchmarks (reference and optimised kernels)",
                 "stratonix bench [description.yaml] [--spacing 50,25] [--mode both|reference|optimised] "
                 "[--voxel 200,25] [--threads N] [-o <folder>] | stratonix bench --write-model <file> [--spacing m]",
                 cli::cmd_bench}},
      {"suite", {"Run the model-interface conformance suite on a description's model",
                 "stratonix suite <description.yaml>", cli::cmd_suite}},
      {"verify", {"Verify a manifest.json against the files beside it", "stratonix verify <manifest.json>",
                  cli::cmd_verify}},
      {"schema", {"Print the JSON Schema of model descriptions", "stratonix schema", cli::cmd_schema}},
      {"version", {"Print version and build information", "stratonix version", cli::cmd_version}},
  };
  return c;
}

void usage() {
  std::cout << "Stratonix " << build_info().version << " - geomodelling toolkit\n\n"
            << "usage: stratonix <command> [options]\n\ncommands:\n";
  for (auto& [n, c] : commands()) {
    std::string pad(12 - std::min<std::size_t>(11, n.size()), ' ');
    std::cout << "  " << n << pad << c.summary << "\n";
  }
  std::cout << "\nRun 'stratonix <command> --help' for the options of a command.\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> a(argv + 1, argv + argc);
  if (a.empty() || a[0] == "-h" || a[0] == "--help" || a[0] == "help") {
    usage();
    return a.empty() ? 2 : 0;
  }
  auto it = commands().find(a[0]);
  if (it == commands().end()) {
    std::cerr << "stratonix: unknown command '" << a[0] << "'\n";
    usage();
    return 2;
  }
  cli::Args args(std::vector<std::string>(a.begin() + 1, a.end()));
  if (args.flag("help") || args.flag("h")) {
    std::cout << "usage: " << it->second.usage << "\n";
    return 0;
  }
  try {
    register_all_plugins();
    if (args.has("threads")) set_thread_count(static_cast<unsigned>(std::stoul(args.get("threads"))));
    if (args.flag("verbose")) set_log_level(LogLevel::Debug);
    return it->second.run(args);
  } catch (const Error& e) {
    std::cerr << "stratonix " << a[0] << ": " << e.what() << "\n";
    switch (e.code()) {
      case ErrorCode::Validation: return 3;
      case ErrorCode::CheckFailed: return 4;
      case ErrorCode::NotImplemented: return 5;
      default: return 1;
    }
  } catch (const std::exception& e) {
    std::cerr << "stratonix " << a[0] << ": internal error: " << e.what() << "\n";
    return 1;
  }
}
