// Delivery package: README, parameters table, checks report, exports, figures, manifest, zip,
// banned-pattern scan and input sharing status.
#include "stratonix/package/delivery.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <regex>
#include <set>
#include <sstream>

#include "stratonix/checks/checks.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/desc/document.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/figures/canvas.hpp"
#include "stratonix/io/export.hpp"

namespace stratonix {

namespace fs = std::filesystem;

// Defined in pipeline.cpp.
void ensure_new_output_dir(const fs::path& out_dir);
BuildOutcome build_only(const desc::Description& d, const BuildOptions& opt);
std::vector<io::ExportedFile> write_exports(const desc::Description& d, const Model& m, const std::vector<std::string>& ids, const fs::path& folder,
                                            double grid_spacing, const GridProducts* gp = nullptr, repr::FaultedSurfaceCache* surfaces = nullptr);
std::vector<std::pair<std::string, std::string>> grid_rows(const GridProducts& gp);
void write_grid_outputs(const desc::Description& d, const Model& m, const GridProducts& gp, const std::vector<std::pair<std::string, std::string>>& keys,
                        const fs::path& out_dir, bool figures);

std::string evidence_in_plain_words(desc::Evidence e) {
  switch (e) {
    case desc::Evidence::Measured: return "measured directly (for example in a well)";
    case desc::Evidence::Interpreted: return "interpreted from data by a geoscientist";
    case desc::Evidence::Calibrated: return "adjusted so the model matches observations";
    case desc::Evidence::Upscaled: return "derived from finer-scale data";
    case desc::Evidence::LiteratureDerived: return "taken from published literature";
    case desc::Evidence::Assumed: return "assumed, with no direct data behind it";
    case desc::Evidence::UserInstruction: return "set on the instruction of the model owner (a design choice, not a measurement)";
  }
  return "unknown";
}

// ---- configuration ---------------------------------------------------------------------------
std::vector<std::string> PackageConfig::builtin_patterns() {
  return {R"(\b[DIOS]-[0-9]{3}\b)",              // register entries: decisions, issues, outputs, scripts
          R"(\bR-[0-9]+\.[0-9]+\.[0-9]+\b)",      // register release entries
          R"(\bWBS-[A-Z]+)",                      // work-breakdown codes
          R"(\bSTX-(TN|RI)-[0-9]+)",              // internal note identifiers
          R"([A-Za-z]:\\[A-Za-z0-9_])",           // absolute Windows paths
          R"(/(home|Users|mnt|tmp)/[A-Za-z0-9_])"};  // absolute Unix paths
}

std::vector<std::pair<std::string, std::string>> PackageConfig::builtin_translations() {
  return {{"Petrel import checklist (STX-TN-02)", "Petrel import checklist supplied with the package documentation"},
          {"(STX-TN-02)", "(the Petrel import checklist)"}};
}

PackageConfig PackageConfig::load(const fs::path& file) { return parse(read_file(file), file.string()); }

PackageConfig PackageConfig::parse(const std::string& text, const std::string& origin) {
  PackageConfig c;
  auto n = desc::parse_document(text, origin);
  if (!n.is_map()) throw Error(ErrorCode::Parse, "package config must be a mapping: " + origin);
  auto strs = [&](const char* k, std::vector<std::string>& out) {
    if (auto* v = n.get(k)) {
      out.clear();
      for (auto& x : v->seq) out.push_back(x.str);
    }
  };
  strs("formats", c.formats);
  strs("assumptions", c.assumptions);
  strs("limitations", c.limitations);
  strs("banned_patterns", c.banned_patterns);
  auto flag = [](const desc::Node& v) { return v.type == desc::Node::Type::Bool ? v.boolean : v.str == "true"; };
  if (auto* v = n.get("figures")) c.figures = flag(*v);
  if (auto* v = n.get("zip")) c.zip = flag(*v);
  if (auto* v = n.get("default_patterns")) c.default_patterns = flag(*v);
  if (auto* v = n.get("terms_of_use")) c.terms_of_use = v->str;
  if (auto* v = n.get("contact")) c.contact = v->str;
  if (auto* v = n.get("translations"))
    for (auto& kv : v->map) c.translations.push_back({kv.first, kv.second.str});
  for (auto& p : c.banned_patterns) {
    try {
      std::regex r(p);
    } catch (const std::regex_error&) {
      throw Error(ErrorCode::Validation, "package config: invalid regular expression in banned_patterns: " + p);
    }
  }
  return c;
}

// ---- scan --------------------------------------------------------------------------------------
namespace {
bool is_text(const fs::path& p, const std::string& data) {
  static const std::set<std::string> binary = {".png", ".tif", ".tiff", ".pptx", ".zip", ".pdf"};
  if (binary.count(p.extension().string())) return false;
  return data.find('\0') == std::string::npos;
}
}  // namespace

std::vector<ScanHit> scan_folder(const fs::path& root, const std::vector<std::string>& patterns) {
  std::vector<std::regex> rx;
  for (auto& p : patterns) rx.emplace_back(p);
  std::vector<ScanHit> hits;
  std::vector<fs::path> files;
  for (auto& e : fs::recursive_directory_iterator(root)) files.push_back(e.path());
  std::sort(files.begin(), files.end());
  for (auto& f : files) {
    std::string rel = fs::relative(f, root).generic_string();
    for (std::size_t k = 0; k < rx.size(); ++k)
      if (std::regex_search(rel, rx[k])) hits.push_back({rel, patterns[k], "(file name)", 0});
    if (!fs::is_regular_file(f)) continue;
    std::string data = read_file(f);
    if (!is_text(f, data)) continue;
    std::size_t line = 1, start = 0;
    for (std::size_t i = 0; i <= data.size(); ++i) {
      if (i == data.size() || data[i] == '\n') {
        std::string l = data.substr(start, i - start);
        for (std::size_t k = 0; k < rx.size(); ++k) {
          std::smatch m;
          if (std::regex_search(l, m, rx[k])) {
            auto p = static_cast<std::size_t>(m.position(0));
            std::size_t a = p > 30 ? p - 30 : 0;
            hits.push_back({rel, patterns[k], l.substr(a, 80), line});
          }
        }
        ++line;
        start = i + 1;
      }
    }
  }
  return hits;
}

// ---- package -----------------------------------------------------------------------------------
namespace {
std::string csv_field(const std::string& s) {
  if (s.find_first_of(",\"\n") == std::string::npos) return s;
  std::string o = "\"";
  for (char c : s) o += (c == '"') ? std::string("\"\"") : std::string(1, c);
  return o + "\"";
}

std::string wrap(const std::string& text, std::size_t width, const std::string& indent) {
  std::string out, line;
  std::istringstream in(text);
  std::string w;
  while (in >> w) {
    if (!line.empty() && line.size() + 1 + w.size() > width) {
      out += indent + line + "\n";
      line.clear();
    }
    line += (line.empty() ? "" : " ") + w;
  }
  if (!line.empty()) out += indent + line + "\n";
  return out;
}

std::string parameters_csv(const desc::Description& d) {
  std::string o = "parameter,value,unit,evidence_class,evidence_in_plain_words,source\n";
  for (auto* p : d.all_params()) {
    std::string val = p->text;
    auto bar = val.find('|');
    if (bar != std::string::npos) val = trim(val.substr(0, bar));
    o += csv_field(p->path) + "," + csv_field(val) + "," + csv_field(p->unit) + "," + desc::evidence_name(p->evidence) + "," +
         csv_field(evidence_in_plain_words(p->evidence)) + "," + csv_field(p->source) + "\n";
  }
  return o;
}

std::string readme(const desc::Description& d, const Model& m, const PackageConfig& cfg, const CheckReport& rep,
                   const std::vector<io::ExportedFile>& files, bool figures) {
  const auto& inf = m.info();
  std::string s;
  bool final_terms = !cfg.terms_of_use.empty();
  s += "Model delivery package\n======================\n\n";
  s += "Model:               " + inf.name + " version " + inf.version + "\n";
  if (!inf.title.empty()) s += "Title:               " + inf.title + "\n";
  s += "Status:              " + std::string(basis_label(inf.basis)) +
       (inf.basis == EvidenceBasis::Design ? " - a design model built from stated parameters, not an interpretation of data"
                                           : " - built from interpreted data") + "\n";
  s += "Package status:      " + std::string(final_terms ? "FINAL" : "DRAFT - terms of use not yet supplied (see below)") + "\n";
  s += "Coordinate system:   " + inf.crs.label() + (inf.crs.is_local() ? " (local frame, no geodetic reference)" : "") + "\n";
  s += "Vertical datum:      " + inf.crs.vertical_datum + "\n";
  s += "Z convention:        elevation, negative below the vertical datum, in every exported file; CSV files also give depth (positive down)\n";
  s += "Units:               metres\n";
  s += "Grid spacing:        " + fmt_shortest(d.grid_spacing.si()) + " m (exported grids and checks)\n";
  s += "Written by:          Stratonix " + build_info().version + "\n\n";

  s += "Horizons (top to base)\n----------------------\n";
  for (std::size_t h = 0; h < inf.horizons.size(); ++h) s += "  " + std::to_string(h + 1) + "  " + inf.horizons[h].name + "\n";
  s += "\nUnits\n-----\n";
  for (auto& u : inf.units) s += "  " + u.name + "\n";
  if (!inf.faults.empty()) {
    s += "\nFaults (oldest first)\n---------------------\n";
    for (std::size_t f = 0; f < inf.faults.size(); ++f) s += "  " + std::to_string(f + 1) + "  " + inf.faults[f].name + "\n";
  }

  s += "\nEvidence behind the parameters\n------------------------------\n";
  std::map<desc::Evidence, int> count;
  for (auto* p : d.all_params()) ++count[p->evidence];
  s += wrap("Every numeric parameter of the model carries a value, a unit, an evidence class and a source; the full list is in "
            "parameters.csv. The classes used here mean:", 96, "");
  for (auto& [e, n] : count)
    s += "  - " + std::string(desc::evidence_name(e)) + " (" + std::to_string(n) + (n == 1 ? " parameter" : " parameters") +
         "): " + evidence_in_plain_words(e) + "\n";

  if (!inf.provenance.empty()) {
    s += "\nData the model honours\n----------------------\n";
    for (auto& p : inf.provenance) s += wrap(p, 94, "    ").replace(0, 4, "  - ");
    s += "  Misfit at every data point, with the evidence class of its input: exports/data/data_misfit.csv\n";
  }
  s += "\nAssumptions\n-----------\n";
  std::vector<std::string> as = cfg.assumptions;
  if (count.count(desc::Evidence::Assumed))
    for (auto* p : d.all_params())
      if (p->evidence == desc::Evidence::Assumed)
        as.push_back("Assumed: " + p->path + " = " + trim(p->text.substr(0, p->text.find('|'))) + (p->source.empty() ? "" : " (" + p->source + ")") + ".");
  if (inf.basis == EvidenceBasis::Design) as.push_back("The model is a design case: its geometry follows the stated parameters and is not fitted to data beyond the control wells.");
  if (as.empty()) as.push_back("None recorded beyond the evidence classes above.");
  for (auto& a : as) s += wrap(a, 94, "    ").replace(0, 4, "  - ");

  s += "\nKnown limitations\n-----------------\n";
  std::vector<std::string> lim = cfg.limitations;
  lim.push_back("Faults are planar; throw follows the stated profiles. The model does not represent fault-zone width or drag.");
  lim.push_back("Petrel import behaviour of the exported formats is to be confirmed with the Petrel import checklist.");
  lim.push_back("Raster figures (PNG, TIFF, PPTX) depend on the fonts of the machine that drew them; the SVG figures carry the same content as text.");
  for (auto& a : lim) s += wrap(a, 94, "    ").replace(0, 4, "  - ");

  s += "\nChecks\n------\n";
  s += wrap(std::to_string(rep.results.size()) + " checks, " + std::to_string(rep.failures()) + " failed (" + (rep.passed() ? "PASS" : "FAIL") +
                "). Every check with its value, tolerance and result is in checks_report.txt and checks_report.csv.", 96, "  ");

  s += "\nContents\n--------\n";
  s += "  README.txt              this file\n  parameters.csv          every parameter: value, unit, evidence class (and in plain words), source\n";
  s += "  model_description.yaml  the model description the package was built from\n";
  s += "  checks_report.txt/.csv  every check, its value, tolerance and pass or fail\n  summary_table.csv       key numbers\n";
  s += "  exports/                model files (" + std::to_string(files.size()) + " files; exports/README.txt lists each with its import option)\n";
  if (figures) s += "  figures/                sections, maps, 3D views and fence diagram in SVG, PNG, TIFF and PPTX\n";
  s += "  manifest.json           SHA-256 of every file; verify with: stratonix verify manifest.json\n";

  s += "\nTerms of use\n------------\n";
  if (final_terms) s += wrap(cfg.terms_of_use, 96, "  ");
  else s += wrap("Not supplied. This package is a DRAFT and must not be shared outside the model owner's organisation until terms of use "
                 "approved by the model owner are added (package configuration: terms_of_use).", 96, "  ");
  s += "\nContact\n-------\n";
  std::string contact = !cfg.contact.empty() ? cfg.contact : d.model.contact;
  s += "  " + (contact.empty() ? std::string("Not supplied (package configuration: contact, or model.contact in the description).") : contact) + "\n";
  return s;
}

std::string apply_translations(std::string s, const std::vector<std::pair<std::string, std::string>>& tr) {
  for (auto& [a, b] : tr) {
    if (a.empty()) continue;
    for (std::size_t p = 0; (p = s.find(a, p)) != std::string::npos; p += b.size()) s.replace(p, a.size(), b);
  }
  return s;
}
}  // namespace

PackageOutcome run_package(const fs::path& description, const fs::path& out_dir, const PackageConfig& cfg, const BuildOptions& opt) {
  PackageOutcome po;
  auto d = desc::load_description(description);
  // 1. Inputs and their sharing status.
  bool unknown = false;
  for (auto& in : d.inputs) {
    po.input_status.push_back(in.id + ": " + in.sharing + " (" + in.path + ", " + in.role + ")");
    if (in.sharing == "unknown") unknown = true;
  }
  if (unknown)
    throw Error(ErrorCode::Validation, "an input has sharing status 'unknown'; record its status (internal, shareable or restricted) before packaging",
                "inputs");
  ensure_new_output_dir(out_dir);
  // 2. Build and check.
  BuildOutcome o = build_only(d, opt);
  CheckContext ctx(*o.model, d, make_well_paths(d), d.grid_spacing.si());
  CheckReport rep = run_checks(ctx);
  po.checks_passed = rep.passed();
  if (!po.checks_passed) {
    write_file(out_dir / "checks_report.txt", rep.text());
    throw Error(ErrorCode::CheckFailed, "checks failed; no package written (see checks_report.txt in " + out_dir.string() + ")");
  }
  fs::create_directories(out_dir);
  write_file(out_dir / "checks_report.txt", rep.text());
  write_file(out_dir / "checks_report.csv", rep.csv());
  auto keys = key_numbers(d, *o.model, rep, ctx);
  const GridProducts* gp = ctx.grid_products();
  {
    auto rows = keys;
    if (gp)
      for (auto& r : grid_rows(*gp)) rows.push_back(r);
    write_file(out_dir / "summary_table.csv", summary_csv(rows));
  }
  write_file(out_dir / "parameters.csv", parameters_csv(d));
  write_file(out_dir / "model_description.yaml", read_file(description));
  auto ids = resolve_export_ids(d, cfg.formats);
  auto files = write_exports(d, *o.model, ids, out_dir / "exports", 0, gp, ctx.surfaces());
  if (cfg.figures) write_figures(d, *o.model, keys, out_dir / "figures", ctx.surfaces());
  if (gp) write_grid_outputs(d, *o.model, *gp, keys, out_dir, cfg.figures);
  write_file(out_dir / "README.txt", readme(d, *o.model, cfg, rep, files, cfg.figures));
  // 3. Translate internal labels at the boundary, in every text file.
  auto tr = PackageConfig::builtin_translations();
  tr.insert(tr.begin(), cfg.translations.begin(), cfg.translations.end());
  std::vector<fs::path> all;
  for (auto& e : fs::recursive_directory_iterator(out_dir))
    if (e.is_regular_file()) all.push_back(e.path());
  std::sort(all.begin(), all.end());
  for (auto& f : all) {
    std::string data = read_file(f);
    if (!is_text(f, data)) continue;
    std::string t = apply_translations(data, tr);
    if (t != data) write_file(f, t);
  }
  // 4. Scan.
  auto patterns = cfg.banned_patterns;
  if (cfg.default_patterns) {
    auto b = PackageConfig::builtin_patterns();
    patterns.insert(patterns.end(), b.begin(), b.end());
  }
  po.hits = scan_folder(out_dir, patterns);
  {  // the folder name travels with the package (and is the top folder of the zip)
    std::string top = out_dir.filename().string();
    for (auto& p : patterns)
      if (std::regex_search(top, std::regex(p))) po.hits.push_back({top + "/", p, "(package folder name)", 0});
  }
  po.scan_passed = po.hits.empty();
  if (!po.scan_passed) {
    std::string r = "Identifier scan FAILED: " + std::to_string(po.hits.size()) + " occurrence(s) of banned patterns. Nothing may be shared.\n\n";
    for (auto& h : po.hits) r += h.file + ":" + std::to_string(h.line) + "  pattern " + h.pattern + "  ..." + h.excerpt + "...\n";
    fs::path rejected = out_dir.string() + "_REJECTED_SCAN";
    if (fs::exists(rejected)) throw Error(ErrorCode::Io, "folder exists: " + rejected.string());
    fs::rename(out_dir, rejected);
    write_file(rejected / "SCAN_FAILED.txt", r);
    po.folder = rejected;
    return po;
  }
  // 5. Manifest (portable: no compiler or platform, so it differs between machines only by its timestamp).
  Manifest man = Manifest::from_description(d, "package");
  man.portable = true;
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  po.folder = out_dir;
  // 6. Zip (entries sorted, no timestamps).
  if (cfg.zip) {
    fig::ZipWriter z;
    std::vector<fs::path> files2;
    for (auto& e : fs::recursive_directory_iterator(out_dir))
      if (e.is_regular_file()) files2.push_back(e.path());
    std::sort(files2.begin(), files2.end(), [&](const fs::path& a, const fs::path& b) {
      return fs::relative(a, out_dir).generic_string() < fs::relative(b, out_dir).generic_string();
    });
    std::string top = out_dir.filename().string();
    for (auto& f : files2) z.add(top + "/" + fs::relative(f, out_dir).generic_string(), read_file(f));
    po.zip = out_dir.string() + ".zip";
    write_file(po.zip, z.finish());
  }
  return po;
}

}  // namespace stratonix
