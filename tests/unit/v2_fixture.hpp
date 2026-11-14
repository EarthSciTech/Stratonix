// Shared fixture for the stage 2A analytical tests (STX-TN-05): descriptions as text, build,
// horizon depth queries and build-log values.
#pragma once
#include <cmath>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "../stx_test.hpp"
#include "stratonix/checks/checks.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/parametric.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

namespace v2fx {
using namespace stratonix;

inline const char* kP = " | user_instruction | test";
inline std::string q(const std::string& v) { return "\"" + v + kP + "\""; }
inline const std::string kFree = "{rule: free}";
inline std::string flat_top(const std::string& taper) { return "{type: flat_topped, taper: " + q(taper + " fraction") + "}"; }

/// Description text: horizons at the given depths (H1, H2, ...), control well W1 at [1000, 1000],
/// domain 10 x 10 km x 3 km, 50 m grid. `structure` defaults to flat layers.
inline std::string model_yaml(const std::string& faults, const std::string& structure = "",
                              const std::vector<double>& depths = {1000, 1200, 1400, 1600},
                              const std::string& extra = "") {
  std::string s = "stratonix_schema: 1\nmodel: {name: t, version: \"1.0.0\", status: DESIGN}\n"
                  "crs: {name: local, vertical_datum: datum, length_unit: m}\n"
                  "domain:\n  origin: " + q("[0, 0] m") + "\n  size: " + q("[10000, 10000] m") +
                  "\n  depth_range: " + q("[0, 3000] m") + "\n  grid_spacing: " + q("50 m") + "\n";
  s += "wells:\n  - {name: W1, role: control, head: " + q("[1000, 1000] m") + ", tops: {";
  for (std::size_t h = 0; h < depths.size(); ++h)
    s += (h ? ", H" : "H") + std::to_string(h + 1) + ": " + q(fmt_fixed(depths[h], 1) + " m");
  s += "}}\nstratigraphy:\n  reference: {well: W1}\n  horizons: [";
  for (std::size_t h = 0; h < depths.size(); ++h) s += (h ? ", {name: H" : "{name: H") + std::to_string(h + 1) + "}";
  s += "]\n  units: [";
  for (std::size_t h = 1; h < depths.size(); ++h)
    s += (h > 1 ? ", {name: U" : "{name: U") + std::to_string(h) + ", top: H" + std::to_string(h) + ", base: H" +
         std::to_string(h + 1) + "}";
  s += "]\n";
  s += structure.empty() ? "structure:\n  regional: {type: flat}\n" : structure;
  if (!faults.empty()) s += "faults:\n" + faults + "frame: {fixed: {well: W1}}\n";
  s += extra;
  return s;
}

/// Validate a description; returns the report (and the description when valid).
inline desc::ValidationReport validate(const std::string& text, desc::Description* d = nullptr) {
  register_all_plugins();
  desc::Description tmp;
  return desc::validate_description_text(text, stxtest::tmp_dir(), d ? d : &tmp);
}

/// Build a description; prints the validation report and throws when it is refused, so an
/// expected-to-fail test states why it fails.
inline ModelPtr build_text(const std::string& text, desc::Description* dout = nullptr, BuildLog* log = nullptr) {
  desc::Description d;
  auto rep = validate(text, &d);
  if (!rep.ok()) {
    std::cerr << "  refused by validation (expected until the stage 2A prompt is done):\n" << rep.text();
    throw std::runtime_error("description refused");
  }
  BuildOptions o;
  o.builder = "parametric";
  auto r = build_only(d, o);
  if (dout) *dout = d;
  if (log) *log = r.log;
  return r.model;
}

/// Depth (positive down) of horizon h in fault block `block` at (x, y).
inline double depth_in_block(const Model& m, int h, BlockCode block, double x, double y) {
  auto g = repr::horizon_grid_block(m, h, block, GridSpec{x, y, 1, 1, 1, 1});
  return -g.v[0];
}
/// Depth of the present-day (faulted) horizon h at (x, y); NaN in a gap.
inline double depth_faulted(const Model& m, int h, double x, double y) {
  auto f = repr::horizon_grid_faulted(m, h, GridSpec{x, y, 1, 1, 1, 1});
  return -f.z.at(0, 0);
}

/// Value of a build-log key, or NaN (and a failed check) when it is missing.
inline double log_value(const BuildLog& log, const std::string& key) {
  for (auto& [k, v] : log.values)
    if (k == key) return std::stod(v);
  std::cerr << "  build log has no value '" << key << "'\n";
  return std::nan("");
}

inline CheckReport checks_for(const Model& m, const desc::Description& d) {
  CheckContext ctx(m, d, make_well_paths(d), 100);
  return run_checks(ctx);
}
inline const CheckResult* find_check(const CheckReport& r, const std::string& id, const std::string& subject = "") {
  for (auto& c : r.results)
    if (c.id == id && (subject.empty() || c.subject.find(subject) != std::string::npos)) return &c;
  return nullptr;
}

/// Planar fault block, with optional extra throw keys (e.g. "sense: reverse") and extra lines.
inline std::string planar_fault(const std::string& name, const std::string& age, const std::string& point,
                                const std::string& depth, const std::string& strike, const std::string& dip,
                                const std::string& length, const std::string& throw_keys, const std::string& mode,
                                const std::string& hw = "{type: rigid}", const std::string& extra = "",
                                const std::string& pstrike = "") {
  return "  - name: " + name + "\n    age: " + q(age + " ordinal") + "\n    geometry: {type: planar, reference_point: " +
         q(point) + ", reference_depth: " + q(depth) + ", strike: " + q(strike) + ", dip: " + q(dip) +
         "}\n    extent: {upper_tip: " + kFree + ", lower_tip: " + kFree + ", length: " + q(length) +
         "}\n    throw: {" + throw_keys + ", profile_dip: " + flat_top("0.25") +
         ", profile_strike: " + (pstrike.empty() ? flat_top("0.01") : pstrike) + "}\n    slip_mode: " + mode +
         "\n    hanging_wall: " + hw + "\n" + extra;
}

/// Bisection for a root of f on [a, b] (f(a) and f(b) of opposite sign).
inline double bisect(const std::function<double(double)>& f, double a, double b) {
  double fa = f(a);
  for (int i = 0; i < 200; ++i) {
    double m = 0.5 * (a + b), fm = f(m);
    if ((fm < 0) == (fa < 0)) { a = m; fa = fm; } else b = m;
  }
  return 0.5 * (a + b);
}

}  // namespace v2fx
