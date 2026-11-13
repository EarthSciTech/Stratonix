// Export context and the Petrel export set (exporter plugins).
#include <algorithm>
#include <cmath>
#include <set>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"

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
    for (std::size_t f = 0; f < inf.faults.size(); ++f) {
      TSurfInfo ti{petrel_name(inf.faults[f].name), "Elevation", 0.8, 0.2, 0.2};
      fs::path p = fs::path("faults") / (petrel_name(inf.faults[f].name) + ".ts");
      const auto& fsu = c.fault_surface(static_cast<int>(f));
      put(folder, p, write_tsurf(fsu.mesh, ti, c.decimals));
      out.push_back({rel(p), "GOCAD TSurf", "GOCAD TSurf (ASCII)",
                     inf.faults[f].name + " triangulated plane clipped to the tip line (" + std::to_string(fsu.mesh.t.size()) + " triangles)"});
    }
    for (std::size_t h = 0; h < inf.horizons.size(); ++h) {
      TSurfInfo ti{petrel_name(inf.horizons[h].name), "Elevation", 0.6, 0.7, 0.9};
      fs::path p = fs::path("horizons") / (petrel_name(inf.horizons[h].name) + "_faulted.ts");
      TriMesh mesh = repr::grid_to_mesh(c.faulted(static_cast<int>(h)).z);
      put(folder, p, write_tsurf(mesh, ti, c.decimals));
      out.push_back({rel(p), "GOCAD TSurf", "GOCAD TSurf (ASCII)",
                     inf.horizons[h].name + " faulted surface as triangles (no triangles across gaps)"});
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

void register_exporters() {
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
  write_file(folder / "README.txt", export_readme(c, all, evidence_note));
  return all;
}

}  // namespace stratonix::io
