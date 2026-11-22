// Export context and the Petrel export set (exporter plugins).
#include <algorithm>
#include <map>
#include <cmath>
#include <set>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/engines/structure.hpp"
#include "stratonix/io/io.hpp"

// Plugin anchors called below (C linkage, defined by STRATONIX_REGISTER in each plugin file).
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_irap_ascii);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_zmap);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_fault_polygons);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_fault_sticks);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_petrel_points);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_tsurf);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_well_heads);
STRATONIX_DECLARE_PLUGIN_ANCHOR(exp_well_tops);

namespace stratonix::io {

namespace fs = std::filesystem;

// ---- context ---------------------------------------------------------------------------------
ExportContext::ExportContext(const Model& m, std::vector<WellPath> wells, const GridSpec& grid, double fault_spacing, double stick_spacing)
    : m_(m), wells_(std::move(wells)), grid_(grid), fault_spacing_(fault_spacing), stick_spacing_(stick_spacing) {}

std::vector<std::string> ExportContext::comments() const {
  const auto& i = info();
  return {"Model: " + i.name + " v" + i.version + " (" + basis_label(i.basis) + ")",
          "Coordinate system: " + i.crs.label() + "; vertical datum: " + i.crs.vertical_datum,
          "Units: m. Z is elevation (negative below the datum).", "Written by Stratonix " + build_info().version};
}

const repr::FaultedSurface& ExportContext::faulted(int h) {
  if (shared_ && &shared_->model() == &m_) return shared_->get(h, grid_);
  auto it = faulted_.find(h);
  if (it == faulted_.end()) it = faulted_.emplace(h, repr::horizon_grid_faulted(m_, h, grid_)).first;
  return it->second;
}
const repr::FaultSurface& ExportContext::fault_surface(int f) {
  auto it = fsurf_.find(f);
  if (it == fsurf_.end()) it = fsurf_.emplace(f, repr::fault_surface(m_, f, fault_spacing_)).first;
  return it->second;
}
const std::vector<repr::CutoffLine>& ExportContext::cutoffs(int f, int h) {
  auto key = std::make_pair(f, h);
  auto it = cut_.find(key);
  if (it == cut_.end()) it = cut_.emplace(key, repr::cutoff_lines(m_, fault_surface(f), h)).first;
  return it->second;
}
const std::vector<repr::WellTop>& ExportContext::tops(std::size_t w) {
  auto it = tops_.find(w);
  if (it == tops_.end()) it = tops_.emplace(w, repr::well_intersections(m_, wells_.at(w), 2.0)).first;
  return it->second;
}

std::vector<BlockCode> ExportContext::blocks_on(int h) {
  const auto& fsu = faulted(h);
  std::set<BlockCode> s;
  for (std::size_t k = 0; k < fsu.z.v.size(); ++k)
    if (!std::isnan(fsu.z.v[k])) s.insert(fsu.block[k]);
  return {s.begin(), s.end()};
}

Grid2D ExportContext::block_surface(int h, BlockCode b, int dilate) {
  const auto& fsu = faulted(h);
  const auto& g = grid_;
  std::vector<unsigned char> mask(g.size(), 0);
  for (std::size_t k = 0; k < g.size(); ++k) mask[k] = !std::isnan(fsu.z.v[k]) && fsu.block[k] == b;
  for (int it = 0; it < dilate; ++it) {
    auto m2 = mask;
    for (std::size_t j = 0; j < g.ny; ++j)
      for (std::size_t i = 0; i < g.nx; ++i) {
        if (mask[g.index(i, j)]) continue;
        bool nb = (i > 0 && mask[g.index(i - 1, j)]) || (i + 1 < g.nx && mask[g.index(i + 1, j)]) || (j > 0 && mask[g.index(i, j - 1)]) ||
                  (j + 1 < g.ny && mask[g.index(i, j + 1)]);
        if (nb) m2[g.index(i, j)] = 1;
      }
    mask.swap(m2);
  }
  Grid2D full = repr::horizon_grid_block(m_, h, b, g);
  bool can_overlap = false;
  for (auto& fi : info().faults) can_overlap = can_overlap || fi.overlapping;
  if (can_overlap) {
    // Under a reverse fault the faulted surface keeps the upper sheet; the lower block still holds
    // the horizon there, so its own surface covers every node where the model puts it in the block.
    std::vector<std::size_t> idx;
    std::vector<Vec3> p;
    for (std::size_t k = 0; k < g.size(); ++k)
      if (!mask[k] && !std::isnan(full.v[k]) && full.v[k] <= info().domain.z_top && full.v[k] >= info().domain.z_base) {
        idx.push_back(k);
        p.push_back({g.x(k % g.nx), g.y(k / g.nx), full.v[k]});
      }
    EvalResult r;
    m_.evaluate(p, kEvalBlock, r);
    for (std::size_t q = 0; q < idx.size(); ++q)
      if (r.block[q] == b) mask[idx[q]] = 1;
  }
  for (std::size_t k = 0; k < g.size(); ++k)
    if (!mask[k]) full.v[k] = kNaN;
  return full;
}

std::string ExportContext::block_label(BlockCode b) const {
  std::string s;
  for (std::size_t f = 0; f < info().faults.size(); ++f) {
    if (!s.empty()) s += "_";
    s += petrel_name(info().faults[f].name) + (((b >> f) & 1u) ? "-HW" : "-FW");
  }
  return s.empty() ? "single" : s;
}

std::vector<Polyline> ExportContext::fault_polygons(int f, int h) {
  const auto& lines = cutoffs(f, h);
  const repr::CutoffLine* hw = nullptr;
  const repr::CutoffLine* fw = nullptr;
  for (auto& l : lines) (l.hanging_wall ? hw : fw) = &l;
  std::vector<Polyline> out;
  if (!hw || !fw) return out;
  Polyline a, b;
  auto flush = [&] {
    if (a.size() >= 2) {
      Polyline poly = a;
      poly.insert(poly.end(), b.rbegin(), b.rend());
      poly.push_back(poly.front());
      out.push_back(poly);
    }
    a.clear();
    b.clear();
  };
  // Runs of columns where both cut-offs exist; a break of up to four columns (a cut-off missed
  // near a tip or where a younger fault crosses) is bridged so one gap gives one polygon.
  std::size_t missing = 0;
  for (std::size_t i = 0; i < hw->by_column.size(); ++i) {
    const Vec3 &p = hw->by_column[i], &q = fw->by_column[i];
    if (std::isnan(p.x) || std::isnan(q.x)) {
      if (++missing > 4) flush();
      continue;
    }
    missing = 0;
    a.push_back(q);
    b.push_back(p);
  }
  flush();
  return out;
}

// ---- exporters -------------------------------------------------------------------------------
namespace {
std::string rel(const fs::path& p) { return p.generic_string(); }

void put(const fs::path& folder, const fs::path& relpath, const std::string& text) {
  fs::create_directories((folder / relpath).parent_path());
  write_file(folder / relpath, text);
}

class IrapGridExporter final : public Exporter {
 public:
  std::string id() const override { return "irap_ascii"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      std::string hn = petrel_name(inf.horizons[h].name);
      fs::path p = fs::path("horizons") / (hn + "_faulted.irap");
      put(folder, p, write_irap_ascii(c.faulted(static_cast<int>(h)).z, c.decimals));
      out.push_back({rel(p), "IRAP classic ASCII grid", "IRAP classic grid (ASCII)",
                     inf.horizons[h].name + " faulted surface; fault gaps undefined (9999900.0)"});
      // Per-block surfaces only where a fault offsets the horizon (it has fault polygons).
      bool cut = false;
      for (std::size_t f = 0; f < inf.faults.size() && !cut; ++f) cut = !c.fault_polygons(static_cast<int>(f), static_cast<int>(h)).empty();
      if (!cut) continue;
      for (BlockCode b : c.blocks_on(static_cast<int>(h))) {
        fs::path pb = fs::path("horizons") / "per_block" / (hn + "_block" + std::to_string(b) + ".irap");
        put(folder, pb, write_irap_ascii(c.block_surface(static_cast<int>(h), b), c.decimals));
        out.push_back({rel(pb), "IRAP classic ASCII grid", "IRAP classic grid (ASCII)",
                       inf.horizons[h].name + " in block " + std::to_string(b) + " (" + c.block_label(b) + "), extended 2 nodes across the gap"});
      }
    }
    // Stage 2D: present-day mean porosity of each compacting unit (for stage 3B), restored-frame
    // map positions.
    if (auto* st = static_cast<const Structure*>(c.model().capability("structure")); st && st->compaction_mode() != 0)
      for (std::size_t u = 1; u + 1 < inf.units.size(); ++u) {
        if (!st->compaction_law(static_cast<int>(u))) continue;
        Grid2D g(c.grid(), std::nan(""));
        for (std::size_t k = 0; k < g.v.size(); ++k) g.v[k] = st->porosity(static_cast<int>(u), c.grid().x(k % c.grid().nx), c.grid().y(k / c.grid().nx));
        fs::path p = fs::path("porosity") / (petrel_name(inf.units[u].name) + "_porosity.irap");
        put(folder, p, write_irap_ascii(g, c.decimals));
        out.push_back({rel(p), "IRAP classic ASCII grid", "IRAP classic grid (ASCII)",
                       inf.units[u].name + " present-day mean porosity (fraction) from its compaction law; restored-frame positions"});
      }
    // Stage 2F: curvature attributes of the horizons asked for (STX-TN-10 section 5).
    for (int h : c.curvature_horizons) {
      auto cg = repr::curvature(c.model(), h, c.grid(), std::isnan(c.fracture_threshold) ? 1e300 : c.fracture_threshold);
      const std::string hn = petrel_name(inf.horizons[static_cast<std::size_t>(h)].name);
      for (auto [grid, name, what] : {std::tuple<const Grid2D*, const char*, const char*>{&cg.mean, "mean", "mean curvature (1/m)"},
                                      {&cg.gaussian, "gaussian", "Gaussian curvature (1/m^2)"},
                                      {&cg.kmax, "kmax", "maximum principal curvature (1/m)"},
                                      {&cg.kmin, "kmin", "minimum principal curvature (1/m)"},
                                      {&cg.most_pos, "most_positive", "most positive curvature (1/m)"},
                                      {&cg.most_neg, "most_negative", "most negative curvature (1/m)"},
                                      {&cg.dip, "dip_curvature", "dip curvature (1/m)"},
                                      {&cg.strike, "strike_curvature", "strike curvature (1/m)"},
                                      {&cg.trend, "fracture_trend", "fracture-trend azimuth (deg) where the curvature passes the threshold"}}) {
        fs::path p = fs::path("curvature") / (hn + "_" + name + ".irap");
        put(folder, p, write_irap_ascii(*grid, std::max(c.decimals, 9)));
        out.push_back({rel(p), "IRAP classic ASCII grid", "IRAP classic grid (ASCII)",
                       inf.horizons[static_cast<std::size_t>(h)].name + " " + what + "; convex up positive; fault gaps undefined"});
      }
    }
    // Stage 2C: isopachs and stratal termination lines of units that terminate (STX-TN-07 section 2).
    for (std::size_t u = 1; u + 1 < inf.units.size(); ++u) {
      if (inf.units[u].termination.empty()) continue;
      std::string un = petrel_name(inf.units[u].name);
      fs::path p = fs::path("isopachs") / (un + "_isopach.irap");
      put(folder, p, write_irap_ascii(repr::unit_thickness_grid(c.model(), static_cast<int>(u), c.grid()), c.decimals));
      out.push_back({rel(p), "IRAP classic ASCII grid", "IRAP classic grid (ASCII)",
                     inf.units[u].name + " isopach (vertical thickness, m); fault gaps undefined"});
      std::map<std::string, std::vector<Polyline>> by_kind;
      for (auto& t : repr::termination_lines(c.model(), static_cast<int>(u), c.grid())) {
        Polyline l = t.line;
        for (auto& q : l) q.z = 0;
        by_kind[t.kind].push_back(l);
      }
      for (auto& [kind, lines] : by_kind) {
        fs::path pl = fs::path("lines") / (un + "_" + kind + "_line.irap");
        put(folder, pl, write_irap_lines(lines, c.decimals));
        out.push_back({rel(pl), "IRAP classic lines", "IRAP classic lines (ASCII)",
                       inf.units[u].name + " " + kind + " line (zero thickness), map position, z = 0"});
      }
    }
    return out;
  }
};

class ZmapExporter final : public Exporter {
 public:
  std::string id() const override { return "zmap"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      std::string hn = petrel_name(inf.horizons[h].name);
      fs::path p = fs::path("horizons") / (hn + "_faulted.zmap");
      put(folder, p, write_zmap(c.faulted(static_cast<int>(h)).z, hn, c.comments(), c.decimals));
      out.push_back({rel(p), "ZMAP+ grid", "Zmap+ grid (ASCII)", inf.horizons[h].name + " faulted surface; fault gaps null (-99999.0)"});
    }
    return out;
  }
};

class FaultPolygonExporter final : public Exporter {
 public:
  std::string id() const override { return "fault_polygons_irap"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      std::vector<Polyline> all;
      std::string which;
      for (std::size_t f = 0; f < inf.faults.size(); ++f) {
        auto p = c.fault_polygons(static_cast<int>(f), static_cast<int>(h));
        if (!p.empty()) which += (which.empty() ? "" : ", ") + inf.faults[f].name;
        all.insert(all.end(), p.begin(), p.end());
      }
      if (all.empty()) continue;
      fs::path p = fs::path("horizons") / "fault_polygons" / (petrel_name(inf.horizons[h].name) + "_fault_polygons.dat");
      put(folder, p, write_irap_lines(all, c.decimals));
      out.push_back({rel(p), "IRAP classic lines", "Irap classic lines (ASCII), as fault polygons",
                     inf.horizons[h].name + " fault polygons (" + which + "): footwall cut-off then hanging-wall cut-off, closed"});
    }
    return out;
  }
};

class FaultSticksIrapExporter final : public Exporter {
 public:
  std::string id() const override { return "fault_sticks_irap"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      auto sticks = repr::fault_sticks(c.fault_surface(static_cast<int>(f)), c.stick_spacing());
      fs::path p = fs::path("faults") / (petrel_name(inf.faults[f].name) + "_sticks.dat");
      put(folder, p, write_irap_lines(sticks, c.decimals));
      out.push_back({rel(p), "IRAP classic lines", "Irap classic lines (ASCII), as fault sticks",
                     inf.faults[f].name + ": " + std::to_string(sticks.size()) + " sticks at " + fmt_shortest(c.stick_spacing()) + " m"});
    }
    return out;
  }
};

class PetrelPointsExporter final : public Exporter {
 public:
  std::string id() const override { return "petrel_points"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    const auto& inf = c.info();
    if (inf.faults.empty()) return {};
    PointsTable t;
    t.attr = {{"Fault", "String"}, {"Stick", "Int"}};
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      auto sticks = repr::fault_sticks(c.fault_surface(static_cast<int>(f)), c.stick_spacing());
      for (std::size_t s = 0; s < sticks.size(); ++s)
        for (auto& p : sticks[s]) {
          t.pts.push_back(p);
          t.values.push_back({petrel_name(inf.faults[f].name), std::to_string(s + 1)});
        }
    }
    fs::path p = fs::path("faults") / "fault_sticks_points.txt";
    put(folder, p, write_petrel_points(t, c.comments(), c.decimals));
    return {{rel(p), "Petrel points with attributes", "Petrel points with attributes (ASCII)",
             "fault stick vertices of all faults, attributes Fault and Stick"}};
  }
};

class TSurfExporter final : public Exporter {
 public:
  std::string id() const override { return "gocad_tsurf"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    std::vector<ExportedFile> out;
    const auto& inf = c.info();
    const std::string crs = petrel_name(inf.crs.label());
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      const auto& fsu = c.fault_surface(static_cast<int>(f));
      TSurfInfo ti;
      ti.name = petrel_name(inf.faults[f].name);
      ti.crs_name = crs;
      ti.r = 0.8; ti.g = 0.2; ti.b = 0.2;
      FaultEval fe;
      c.model().evaluate_fault(static_cast<int>(f), fsu.mesh.v, fe);
      ti.property = "displacement";
      ti.values = fe.displacement;
      fs::path p = fs::path("faults") / (petrel_name(inf.faults[f].name) + ".ts");
      put(folder, p, write_tsurf(fsu.mesh, ti, c.decimals));
      out.push_back({rel(p), "GOCAD TSurf", "GOCAD TSurf (ASCII)",
                     inf.faults[f].name + " triangulated plane clipped to the tip line (" + std::to_string(fsu.mesh.t.size()) +
                         " triangles); vertex property displacement (m)"});
    }
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      TriMesh mesh = repr::grid_to_mesh(c.faulted(static_cast<int>(h)).z);
      TSurfInfo ti;
      ti.name = petrel_name(inf.horizons[h].name);
      ti.crs_name = crs;
      ti.r = 0.6; ti.g = 0.7; ti.b = 0.9;
      ti.property = "depth";
      for (auto& v : mesh.v) ti.values.push_back(-v.z);
      fs::path p = fs::path("horizons") / (petrel_name(inf.horizons[h].name) + "_faulted.ts");
      put(folder, p, write_tsurf(mesh, ti, c.decimals));
      out.push_back({rel(p), "GOCAD TSurf", "GOCAD TSurf (ASCII)",
                     inf.horizons[h].name + " faulted surface as triangles (no triangles across gaps); vertex property depth (m, positive down)"});
    }
    return out;
  }
};

class WellHeadsExporter final : public Exporter {
 public:
  std::string id() const override { return "well_heads"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    if (c.wells().empty()) return {};
    std::vector<WellHead> hs;
    std::vector<ExportedFile> out;
    for (auto& w : c.wells()) {
      hs.push_back({petrel_name(w.name()), w.x(), w.y(), w.datum_z(), w.total_depth_md()});
      if (!w.vertical()) {
        fs::path p = fs::path("wells") / (petrel_name(w.name()) + "_path.txt");
        put(folder, p, write_well_path(w));
        out.push_back({rel(p), "Well path MD X Y Z", "Well path/deviation (ASCII)", w.name() + " deviation survey"});
      }
    }
    fs::path p = fs::path("wells") / "well_heads.txt";
    put(folder, p, write_well_heads(hs, c.comments()));
    out.insert(out.begin(), ExportedFile{rel(p), "Petrel well heads", "Well heads", std::to_string(hs.size()) + " wells, KB elevation, TD as MD"});
    return out;
  }
};

class WellTopsExporter final : public Exporter {
 public:
  std::string id() const override { return "well_tops"; }
  std::vector<ExportedFile> write(ExportContext& c, const fs::path& folder) const override {
    if (c.wells().empty()) return {};
    std::vector<WellTopRow> rows;
    for (std::size_t w = 0; w < c.wells().size(); ++w)
      for (auto& t : c.tops(w)) rows.push_back({petrel_name(c.wells()[w].name()), petrel_name(t.horizon), t.p.x, t.p.y, t.p.z, t.md});
    fs::path p = fs::path("wells") / "well_tops.txt";
    put(folder, p, write_well_tops(rows, c.comments()));
    return {{rel(p), "Petrel well tops", "Petrel well tops (ASCII)", std::to_string(rows.size()) + " model tops (horizon crossings) on all wells"}};
  }
};
}  // namespace

STRATONIX_REGISTER(exp_irap_ascii, registry::Kind::Exporter, "irap_ascii", "1.0.0", "Horizons: IRAP classic ASCII grids (faulted with gaps; per block)",
                   [] { return std::make_unique<IrapGridExporter>(); })
STRATONIX_REGISTER(exp_zmap, registry::Kind::Exporter, "zmap", "1.0.0", "Horizons: ZMAP+ grids (faulted with gaps)",
                   [] { return std::make_unique<ZmapExporter>(); })
STRATONIX_REGISTER(exp_fault_polygons, registry::Kind::Exporter, "fault_polygons_irap", "1.0.0", "Fault polygons per horizon: IRAP classic lines",
                   [] { return std::make_unique<FaultPolygonExporter>(); })
STRATONIX_REGISTER(exp_fault_sticks, registry::Kind::Exporter, "fault_sticks_irap", "1.0.0", "Fault sticks: IRAP classic lines",
                   [] { return std::make_unique<FaultSticksIrapExporter>(); })
STRATONIX_REGISTER(exp_petrel_points, registry::Kind::Exporter, "petrel_points", "1.0.0", "Fault sticks: Petrel points with attributes",
                   [] { return std::make_unique<PetrelPointsExporter>(); })
STRATONIX_REGISTER(exp_tsurf, registry::Kind::Exporter, "gocad_tsurf", "1.0.0", "Faults and horizons: GOCAD TSurf",
                   [] { return std::make_unique<TSurfExporter>(); })
STRATONIX_REGISTER(exp_well_heads, registry::Kind::Exporter, "well_heads", "1.0.0", "Wells: Petrel well heads (and paths of deviated wells)",
                   [] { return std::make_unique<WellHeadsExporter>(); })
STRATONIX_REGISTER(exp_well_tops, registry::Kind::Exporter, "well_tops", "1.0.0", "Wells: Petrel well tops",
                   [] { return std::make_unique<WellTopsExporter>(); })

void register_neutral_exporters();
void register_grid_exporters();
void register_exporters() {
  register_neutral_exporters();
  register_grid_exporters();
  STRATONIX_PLUGIN_ANCHOR(exp_irap_ascii);
  STRATONIX_PLUGIN_ANCHOR(exp_zmap);
  STRATONIX_PLUGIN_ANCHOR(exp_fault_polygons);
  STRATONIX_PLUGIN_ANCHOR(exp_fault_sticks);
  STRATONIX_PLUGIN_ANCHOR(exp_petrel_points);
  STRATONIX_PLUGIN_ANCHOR(exp_tsurf);
  STRATONIX_PLUGIN_ANCHOR(exp_well_heads);
  STRATONIX_PLUGIN_ANCHOR(exp_well_tops);
}

std::vector<std::string> petrel_set_ids() {
  return {"irap_ascii", "zmap", "fault_polygons_irap", "fault_sticks_irap", "petrel_points", "gocad_tsurf", "well_heads", "well_tops"};
}

std::string export_readme(const ExportContext& c, const std::vector<ExportedFile>& files, const std::string& evidence_note) {
  const auto& i = c.info();
  const auto& g = c.grid();
  std::string s;
  s += "Stratonix export\n================\n\n";
  s += "Model:              " + i.name + " v" + i.version + (i.title.empty() ? "" : " - " + i.title) + "\n";
  s += std::string("Evidence basis:     ") + basis_label(i.basis) + (i.basis == EvidenceBasis::Design ? " (design model, not an interpretation)" : "") + "\n";
  if (!evidence_note.empty()) s += "                    " + evidence_note + "\n";
  s += "Coordinate system:  " + i.crs.label() + "\n";
  s += "Vertical datum:     " + i.crs.vertical_datum + "\n";
  s += "Units:              metres (x, y, z, MD)\n";
  s += "Z convention:       elevation, negative below the datum, in every file (TSurf: ZPOSITIVE Elevation)\n";
  s += "Grids:              origin (" + fmt_fixed(g.x0, 2) + ", " + fmt_fixed(g.y0, 2) + "), spacing " + fmt_shortest(g.dx) + " x " +
       fmt_shortest(g.dy) + " m, " + std::to_string(g.nx) + " columns x " + std::to_string(g.ny) + " rows, no rotation\n";
  s += "Undefined values:   IRAP 9999900.0, ZMAP+ -99999.0 (fault gaps and outside the horizon)\n";
  s += "Line separator:     IRAP lines use '999.000000 999.000000 999.000000' between lines\n";
  s += "Fault sticks:       every " + fmt_shortest(c.stick_spacing()) + " m along strike\n";
  s += "Written by:         Stratonix " + build_info().version + "\n\n";
  if (!i.provenance.empty()) {
    s += "Data the model honours (engine " + i.engine + ")\n----------------------------\n";
    for (auto& p : i.provenance) s += p + "\n";
    s += "Misfit at every data point, with the evidence class of its input: data/data_misfit.csv\n\n";
  }
  s += "Files\n-----\n";
  std::size_t w = 4;
  for (auto& f : files) w = std::max(w, f.path.size());
  for (auto& f : files) {
    s += f.path + std::string(w + 2 - f.path.size(), ' ') + f.object + "\n";
    s += std::string(w + 2, ' ') + "format: " + f.format + "; Petrel import: " + f.petrel_import + "\n";
  }
  if (!i.faults.empty()) {
    s += "\nFault blocks (per-block surfaces are named <horizon>_block<code>)\n----------------------------------------------------------------\n";
    std::set<BlockCode> all;
    for (auto& f : files) {
      auto p = f.path.rfind("_block");
      if (p != std::string::npos) all.insert(static_cast<BlockCode>(std::stoull(f.path.substr(p + 6))));
    }
    for (BlockCode b : all) s += "block" + std::to_string(b) + ": " + c.block_label(b) + "\n";
  }
  s += "\nThe Petrel import options above are the file types to choose in Petrel's import dialog. Their behaviour\n"
       "in Petrel is to be confirmed with the Petrel import checklist (STX-TN-02).\n";
  return s;
}

std::vector<ExportedFile> export_set(ExportContext& c, const std::vector<std::string>& ids, const fs::path& folder,
                                     const std::string& evidence_note) {
  if (fs::exists(folder) && !fs::is_empty(folder))
    throw Error(ErrorCode::Io, "export folder exists and is not empty (issued outputs are never overwritten): " + folder.string());
  fs::create_directories(folder);
  auto& reg = registry::Registry::global();
  std::vector<ExportedFile> all;
  for (auto& id : ids) {
    auto e = reg.create<Exporter>(registry::Kind::Exporter, id);
    auto files = e->write(c, folder);
    all.insert(all.end(), files.begin(), files.end());
  }
  if (auto* mis = static_cast<const std::vector<DataMisfit>*>(c.model().capability("data_misfit"))) {
    bool extra = std::any_of(mis->begin(), mis->end(), [](const DataMisfit& r) { return !std::isnan(r.sigma) || !std::isnan(r.loo); });
    auto nm = [](double v, int dec) { return std::isnan(v) ? std::string() : fmt_fixed(v, dec); };
    std::string t = std::string("input,evidence,horizon,kind,x,y,z,misfit,misfit_unit") + (extra ? ",sigma,normalised_misfit,leave_one_out_m" : "") + "\n";
    for (auto& r : *mis) {
      t += r.input + "," + r.evidence + "," + r.label + "," + r.kind + "," + fmt_fixed(r.p.x, 2) + "," + fmt_fixed(r.p.y, 2) + "," +
           fmt_fixed(r.p.z, 2) + "," + fmt_fixed(r.misfit, 4) + "," + (r.kind == "orientation" ? "deg" : "m");
      if (extra) t += "," + nm(r.sigma, 4) + "," + nm(r.normalised, 4) + "," + nm(r.loo, 4);
      t += "\n";
    }
    fs::create_directories(folder / "data");
    write_file(folder / "data/data_misfit.csv", t);
    all.push_back({"data/data_misfit.csv", "CSV", "-", "misfit of the model at every data point, with the evidence class of the input"});
  }
  write_file(folder / "README.txt", export_readme(c, all, evidence_note));
  return all;
}

}  // namespace stratonix::io
