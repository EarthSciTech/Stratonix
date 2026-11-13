// Importer plugins for the data v2 needs: points, fault sticks, triangulated surfaces, well tops.
// (Grid importers irap_ascii, zmap and xyz are in grids.cpp.)
#include <map>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"

namespace stratonix::io {

namespace {
class XyzPointsImporter final : public Importer {
 public:
  std::string data_kind() const override { return "points"; }
  ImportedData read(const std::filesystem::path& f) const override {
    ImportedData d;
    d.kind = "points";
    d.points = read_xyz_points(read_file(f));
    return d;
  }
};

class IrapLinesImporter final : public Importer {
 public:
  std::string data_kind() const override { return "lines"; }
  ImportedData read(const std::filesystem::path& f) const override {
    ImportedData d;
    d.kind = "lines";
    d.z_is = "elevation";
    int id = 0;
    for (auto& l : read_irap_lines(read_file(f))) d.lines.push_back({f.stem().string(), ++id, l});
    return d;
  }
};

/// Petrel points with attributes. When the file has a string attribute "Fault" and an integer
/// attribute "Stick" (as Stratonix writes fault sticks), points are grouped into stick lines.
class PetrelPointsImporter final : public Importer {
 public:
  std::string data_kind() const override { return "points"; }
  ImportedData read(const std::filesystem::path& f) const override {
    auto t = read_petrel_points(read_file(f));
    ImportedData d;
    d.z_is = "elevation";
    int fa = -1, sa = -1;
    for (std::size_t a = 0; a < t.attr.size(); ++a) {
      if (t.attr[a].first == "Fault") fa = static_cast<int>(a);
      if (t.attr[a].first == "Stick") sa = static_cast<int>(a);
    }
    if (fa < 0 || sa < 0) {
      d.kind = "points";
      d.points = t.pts;
      return d;
    }
    d.kind = "lines";
    std::map<std::pair<std::string, int>, std::size_t> at;
    for (std::size_t i = 0; i < t.pts.size(); ++i) {
      auto key = std::make_pair(t.values[i][static_cast<std::size_t>(fa)], std::stoi(t.values[i][static_cast<std::size_t>(sa)]));
      auto it = at.find(key);
      if (it == at.end()) {
        it = at.emplace(key, d.lines.size()).first;
        d.lines.push_back({key.first, key.second, {}});
      }
      d.lines[it->second].pts.push_back(t.pts[i]);
    }
    return d;
  }
};

class TSurfImporter final : public Importer {
 public:
  std::string data_kind() const override { return "meshes"; }
  ImportedData read(const std::filesystem::path& f) const override {
    TSurfInfo info;
    ImportedData d;
    d.kind = "meshes";
    auto m = read_tsurf(read_file(f), &info);
    d.z_is = info.z_positive == "Depth" ? "depth" : "elevation";
    d.meshes.push_back({info.name.empty() ? f.stem().string() : info.name, std::move(m)});
    return d;
  }
};

class WellTopsImporter final : public Importer {
 public:
  std::string data_kind() const override { return "tops"; }
  ImportedData read(const std::filesystem::path& f) const override {
    ImportedData d;
    d.kind = "tops";
    d.z_is = "elevation";
    for (auto& r : read_well_tops(read_file(f))) d.tops.push_back({r.well, r.surface, {r.x, r.y, r.z}, r.md});
    return d;
  }
};
}  // namespace

STRATONIX_REGISTER(imp_xyz_points, registry::Kind::Importer, "xyz_points", "1.0.0", "XYZ points (x y z per line)",
                   [] { return std::make_unique<XyzPointsImporter>(); })
STRATONIX_REGISTER(imp_irap_lines, registry::Kind::Importer, "irap_lines", "1.0.0", "IRAP classic lines (fault sticks, polygons)",
                   [] { return std::make_unique<IrapLinesImporter>(); })
STRATONIX_REGISTER(imp_petrel_points, registry::Kind::Importer, "petrel_points", "1.0.0",
                   "Petrel points with attributes (fault sticks when Fault and Stick attributes are present)",
                   [] { return std::make_unique<PetrelPointsImporter>(); })
STRATONIX_REGISTER(imp_tsurf, registry::Kind::Importer, "gocad_tsurf", "1.0.0", "GOCAD TSurf triangulated surface",
                   [] { return std::make_unique<TSurfImporter>(); })
STRATONIX_REGISTER(imp_well_tops, registry::Kind::Importer, "well_tops", "1.0.0", "Petrel well tops (ASCII)",
                   [] { return std::make_unique<WellTopsImporter>(); })

void register_importers() {
  STRATONIX_PLUGIN_ANCHOR(imp_xyz_points);
  STRATONIX_PLUGIN_ANCHOR(imp_irap_lines);
  STRATONIX_PLUGIN_ANCHOR(imp_petrel_points);
  STRATONIX_PLUGIN_ANCHOR(imp_tsurf);
  STRATONIX_PLUGIN_ANCHOR(imp_well_tops);
}

}  // namespace stratonix::io
