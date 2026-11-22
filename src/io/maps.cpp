// Stage 4C (STX-TN-16 section 2.5): geological maps (ESRI shapefile, GeoJSON, OGC GeoPackage),
// outcrop measurements and elevation models (ESRI ASCII grid, GeoTIFF) as model inputs.
#include "stratonix/io/maps.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>

#include "stratonix/core/error.hpp"
#include "stratonix/core/crs.hpp"
#include "stratonix/core/geometry.hpp"
#include "stratonix/desc/document.hpp"

#ifdef STRATONIX_HAVE_SQLITE
#include <sqlite3.h>
#endif

namespace stratonix::io {

namespace {
std::vector<unsigned char> read_bytes(const std::filesystem::path& f) {
  std::ifstream in(f, std::ios::binary);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + f.string());
  return std::vector<unsigned char>((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
std::uint32_t be32(const unsigned char* p) { return (std::uint32_t(p[0]) << 24) | (std::uint32_t(p[1]) << 16) | (std::uint32_t(p[2]) << 8) | p[3]; }
std::uint32_t le32(const unsigned char* p) { return (std::uint32_t(p[3]) << 24) | (std::uint32_t(p[2]) << 16) | (std::uint32_t(p[1]) << 8) | p[0]; }
double led(const unsigned char* p) { std::uint64_t x = 0; for (int i = 7; i >= 0; --i) x = (x << 8) | p[i]; double d; std::memcpy(&d, &x, 8); return d; }
void put_be32(std::vector<unsigned char>& b, std::size_t o, std::uint32_t v) { for (int i = 0; i < 4; ++i) b[o + static_cast<std::size_t>(i)] = static_cast<unsigned char>(v >> (24 - 8 * i)); }
void put_le32(std::vector<unsigned char>& b, std::size_t o, std::uint32_t v) { for (int i = 0; i < 4; ++i) b[o + static_cast<std::size_t>(i)] = static_cast<unsigned char>(v >> (8 * i)); }
void put_led(std::vector<unsigned char>& b, std::size_t o, double d) { std::uint64_t x; std::memcpy(&x, &d, 8); for (int i = 0; i < 8; ++i) b[o + static_cast<std::size_t>(i)] = static_cast<unsigned char>(x >> (8 * i)); }
std::string trim(std::string s) {
  while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.back())) || s.back() == '\0')) s.pop_back();
  std::size_t a = 0;
  while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  return s.substr(a);
}
}  // namespace

// ---- ESRI shapefile (ESRI 1998): .shp geometry, .dbf attributes ---------------------------------
MapLayer read_shapefile(const std::filesystem::path& shp) {
  MapLayer L;
  L.source = shp.string();
  L.name = shp.stem().string();
  const auto b = read_bytes(shp);
  if (b.size() < 100 || be32(b.data()) != 9994) throw Error(ErrorCode::Parse, shp.string() + ": not an ESRI shapefile (file code 9994)");
  std::size_t pos = 100;
  while (pos + 8 <= b.size()) {
    const int num = static_cast<int>(be32(b.data() + pos));
    const std::size_t len = 2 * static_cast<std::size_t>(be32(b.data() + pos + 4));
    pos += 8;
    if (pos + len > b.size()) { L.report.push_back("record " + std::to_string(num) + " is truncated"); break; }
    const unsigned char* p = b.data() + pos;
    MapFeature f;
    f.id = num;
    const int type = static_cast<int>(le32(p));
    const int base = type % 10;  // 1 point, 3 polyline, 5 polygon, 8 multipoint (+10 Z, +20 M)
    f.type = base == 8 ? 1 : base;
    const bool hasz = type > 10 && type < 20;
    if (type == 0) { f.type = 0; }
    else if (base == 1) {
      f.parts.push_back({{led(p + 4), led(p + 12), hasz && len >= 28 ? led(p + 20) : 0.0}});
    } else if (base == 3 || base == 5 || base == 8) {
      const bool multi = base == 8;
      const std::size_t nparts = multi ? 1 : le32(p + 36), npts = le32(p + (multi ? 36 : 40));
      const std::size_t pts0 = multi ? 40 : 44 + 4 * nparts;
      std::vector<std::size_t> starts;
      for (std::size_t k = 0; k < nparts; ++k) starts.push_back(multi ? 0 : le32(p + 44 + 4 * k));
      starts.push_back(npts);
      const std::size_t z0 = pts0 + 16 * npts + 16;
      for (std::size_t k = 0; k < nparts; ++k) {
        std::vector<Vec3> ring;
        for (std::size_t i = starts[k]; i < starts[k + 1]; ++i)
          ring.push_back({led(p + pts0 + 16 * i), led(p + pts0 + 16 * i + 8), hasz && z0 + 8 * npts <= len ? led(p + z0 + 8 * i) : 0.0});
        f.parts.push_back(ring);
      }
    } else {
      L.report.push_back("record " + std::to_string(num) + ": shape type " + std::to_string(type) + " is not read");
    }
    L.features.push_back(std::move(f));
    pos += len;
  }
  // attributes
  std::filesystem::path dbf = shp;
  dbf.replace_extension(".dbf");
  if (!std::filesystem::exists(dbf)) {
    dbf.replace_extension(".DBF");
    if (!std::filesystem::exists(dbf)) { L.report.push_back("no .dbf: features have no attributes"); return L; }
  }
  const auto d = read_bytes(dbf);
  if (d.size() < 32) throw Error(ErrorCode::Parse, dbf.string() + ": not a dBASE file");
  const std::size_t nrec = le32(d.data() + 4), hsize = d[8] | (std::size_t(d[9]) << 8), rsize = d[10] | (std::size_t(d[11]) << 8);
  struct Field { std::string name; std::size_t len; };
  std::vector<Field> fields;
  for (std::size_t o = 32; o + 32 <= hsize && d[o] != 0x0D; o += 32)
    fields.push_back({trim(std::string(reinterpret_cast<const char*>(d.data() + o), 11)), d[o + 16]});
  for (std::size_t r = 0; r < nrec && r < L.features.size(); ++r) {
    std::size_t o = hsize + r * rsize + 1;
    for (auto& fl : fields) {
      if (o + fl.len > d.size()) break;
      L.features[r].attributes[fl.name] = trim(std::string(reinterpret_cast<const char*>(d.data() + o), fl.len));
      o += fl.len;
    }
  }
  if (nrec != L.features.size()) L.report.push_back(std::to_string(nrec) + " attribute records for " + std::to_string(L.features.size()) + " shapes");
  std::filesystem::path prj = shp;
  prj.replace_extension(".prj");
  if (std::filesystem::exists(prj)) {
    const auto t = read_bytes(prj);
    L.crs_wkt.assign(t.begin(), t.end());
  }
  return L;
}

void write_shapefile(const std::filesystem::path& shp, const MapLayer& L) {
  int type = 0;
  for (auto& f : L.features) if (f.type) { type = f.type; break; }
  if (type != 1 && type != 3 && type != 5) throw Error(ErrorCode::InvalidArgument, "write_shapefile: points, polylines or polygons of one type");
  std::vector<unsigned char> body, index;
  double bx0 = 1e300, by0 = 1e300, bx1 = -1e300, by1 = -1e300;
  for (auto& f : L.features)
    for (auto& r : f.parts)
      for (auto& p : r) { bx0 = std::min(bx0, p.x); by0 = std::min(by0, p.y); bx1 = std::max(bx1, p.x); by1 = std::max(by1, p.y); }
  std::size_t off = 100;
  int num = 0;
  for (auto& f : L.features) {
    std::vector<unsigned char> c;
    if (type == 1) {
      c.assign(20, 0);
      put_le32(c, 0, 1);
      put_led(c, 4, f.parts.empty() || f.parts[0].empty() ? 0.0 : f.parts[0][0].x);
      put_led(c, 12, f.parts.empty() || f.parts[0].empty() ? 0.0 : f.parts[0][0].y);
    } else {
      std::size_t np = 0;
      double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
      for (auto& r : f.parts) { np += r.size(); for (auto& p : r) { x0 = std::min(x0, p.x); y0 = std::min(y0, p.y); x1 = std::max(x1, p.x); y1 = std::max(y1, p.y); } }
      c.assign(44 + 4 * f.parts.size() + 16 * np, 0);
      put_le32(c, 0, static_cast<std::uint32_t>(type));
      put_led(c, 4, x0); put_led(c, 12, y0); put_led(c, 20, x1); put_led(c, 28, y1);
      put_le32(c, 36, static_cast<std::uint32_t>(f.parts.size()));
      put_le32(c, 40, static_cast<std::uint32_t>(np));
      std::size_t k = 0, s = 0;
      for (auto& r : f.parts) { put_le32(c, 44 + 4 * k++, static_cast<std::uint32_t>(s)); s += r.size(); }
      std::size_t q = 44 + 4 * f.parts.size();
      for (auto& r : f.parts) for (auto& p : r) { put_led(c, q, p.x); put_led(c, q + 8, p.y); q += 16; }
    }
    std::vector<unsigned char> h(8);
    put_be32(h, 0, static_cast<std::uint32_t>(++num));
    put_be32(h, 4, static_cast<std::uint32_t>(c.size() / 2));
    std::vector<unsigned char> ix(8);
    put_be32(ix, 0, static_cast<std::uint32_t>(off / 2));
    put_be32(ix, 4, static_cast<std::uint32_t>(c.size() / 2));
    index.insert(index.end(), ix.begin(), ix.end());
    body.insert(body.end(), h.begin(), h.end());
    body.insert(body.end(), c.begin(), c.end());
    off += 8 + c.size();
  }
  auto header = [&](std::size_t total) {
    std::vector<unsigned char> h(100, 0);
    put_be32(h, 0, 9994);
    put_be32(h, 24, static_cast<std::uint32_t>(total / 2));
    put_le32(h, 28, 1000);
    put_le32(h, 32, static_cast<std::uint32_t>(type));
    put_led(h, 36, bx0); put_led(h, 44, by0); put_led(h, 52, bx1); put_led(h, 60, by1);
    return h;
  };
  auto write = [](const std::filesystem::path& p, const std::vector<unsigned char>& a, const std::vector<unsigned char>& b) {
    std::ofstream o(p, std::ios::binary);
    o.write(reinterpret_cast<const char*>(a.data()), static_cast<std::streamsize>(a.size()));
    o.write(reinterpret_cast<const char*>(b.data()), static_cast<std::streamsize>(b.size()));
  };
  write(shp, header(100 + body.size()), body);
  std::filesystem::path shx = shp;
  shx.replace_extension(".shx");
  write(shx, header(100 + index.size()), index);
  // dBASE III: character fields as wide as their longest value
  std::vector<std::string> names;
  for (auto& f : L.features) for (auto& [k, v] : f.attributes) if (std::find(names.begin(), names.end(), k) == names.end()) names.push_back(k);
  std::vector<std::size_t> width(names.size(), 1);
  for (auto& f : L.features)
    for (std::size_t i = 0; i < names.size(); ++i) {
      auto it = f.attributes.find(names[i]);
      if (it != f.attributes.end()) width[i] = std::min<std::size_t>(254, std::max(width[i], it->second.size()));
    }
  std::size_t rsize = 1;
  for (auto w : width) rsize += w;
  const std::size_t hsize = 32 + 32 * names.size() + 1;
  std::vector<unsigned char> d(hsize + rsize * L.features.size() + 1, ' ');
  std::fill(d.begin(), d.begin() + static_cast<long>(hsize), 0);
  d[0] = 0x03; d[1] = 126; d[2] = 11; d[3] = 22;
  put_le32(d, 4, static_cast<std::uint32_t>(L.features.size()));
  d[8] = static_cast<unsigned char>(hsize & 0xFF); d[9] = static_cast<unsigned char>(hsize >> 8);
  d[10] = static_cast<unsigned char>(rsize & 0xFF); d[11] = static_cast<unsigned char>(rsize >> 8);
  for (std::size_t i = 0; i < names.size(); ++i) {
    const std::size_t o = 32 + 32 * i;
    for (std::size_t k = 0; k < 10 && k < names[i].size(); ++k) d[o + k] = static_cast<unsigned char>(names[i][k]);
    d[o + 11] = 'C';
    d[o + 16] = static_cast<unsigned char>(width[i]);
  }
  d[hsize - 1] = 0x0D;
  for (std::size_t r = 0; r < L.features.size(); ++r) {
    std::size_t o = hsize + r * rsize;
    d[o++] = ' ';
    for (std::size_t i = 0; i < names.size(); ++i) {
      auto it = L.features[r].attributes.find(names[i]);
      const std::string v = it == L.features[r].attributes.end() ? std::string() : it->second.substr(0, width[i]);
      for (std::size_t k = 0; k < v.size(); ++k) d[o + k] = static_cast<unsigned char>(v[k]);
      o += width[i];
    }
  }
  d.back() = 0x1A;
  std::filesystem::path dbf = shp;
  dbf.replace_extension(".dbf");
  write(dbf, d, {});
}

// ---- GeoJSON (RFC 7946; the 2008 "crs" member is reported) ---------------------------------------
MapLayer read_geojson(const std::filesystem::path& file) {
  MapLayer L;
  L.source = file.string();
  L.name = file.stem().string();
  const auto bytes = read_bytes(file);
  const desc::Node doc = desc::parse_document(std::string(bytes.begin(), bytes.end()), file.string());
  if (const desc::Node* c = doc.get("crs"))
    if (const desc::Node* pr = c->get("properties"))
      if (const desc::Node* nm = pr->get("name")) L.crs_name = nm->str;
  const desc::Node* fs = doc.get("features");
  if (!fs) throw Error(ErrorCode::Parse, file.string() + ": not a GeoJSON feature collection");
  int id = 0;
  auto pt = [](const desc::Node& n) { return Vec3{n.seq.size() > 0 ? n.seq[0].number : 0.0, n.seq.size() > 1 ? n.seq[1].number : 0.0, n.seq.size() > 2 ? n.seq[2].number : 0.0}; };
  for (auto& f : fs->seq) {
    MapFeature m;
    m.id = ++id;
    if (const desc::Node* pr = f.get("properties"))
      for (auto& kv : pr->map) m.attributes[kv.first] = kv.second.is_null() ? std::string() : kv.second.str;
    const desc::Node* g = f.get("geometry");
    const desc::Node* co = g ? g->get("coordinates") : nullptr;
    const std::string type = g && g->get("type") ? g->get("type")->str : "";
    if (!co) { L.report.push_back("feature " + std::to_string(id) + " has no geometry"); continue; }
    if (type == "Point") { m.type = 1; m.parts.push_back({pt(*co)}); }
    else if (type == "MultiPoint") { m.type = 1; for (auto& p : co->seq) m.parts.push_back({pt(p)}); }
    else if (type == "LineString") { m.type = 3; std::vector<Vec3> l; for (auto& p : co->seq) l.push_back(pt(p)); m.parts.push_back(l); }
    else if (type == "MultiLineString" || type == "Polygon") {
      m.type = type == "Polygon" ? 5 : 3;
      for (auto& r : co->seq) { std::vector<Vec3> l; for (auto& p : r.seq) l.push_back(pt(p)); m.parts.push_back(l); }
    } else if (type == "MultiPolygon") {
      m.type = 5;
      for (auto& poly : co->seq) for (auto& r : poly.seq) { std::vector<Vec3> l; for (auto& p : r.seq) l.push_back(pt(p)); m.parts.push_back(l); }
    } else {
      L.report.push_back("feature " + std::to_string(id) + ": geometry " + type + " is not read");
      continue;
    }
    L.features.push_back(std::move(m));
  }
  return L;
}

// ---- GeoPackage (OGC 12-128r19): features through SQLite -----------------------------------------
namespace {
struct Wkb {
  const unsigned char* p;
  std::size_t n, i = 0;
  bool le = true;
  std::uint32_t u32() { std::uint32_t v = le ? le32(p + i) : be32(p + i); i += 4; return v; }
  double dbl() { double d; if (le) d = led(p + i); else { unsigned char t[8]; for (int k = 0; k < 8; ++k) t[k] = p[i + 7 - static_cast<std::size_t>(k)]; d = led(t); } i += 8; return d; }
  void geom(MapFeature& f) {
    le = p[i++] == 1;
    std::uint32_t t = u32();
    bool z = false;
    if (t & 0x80000000u) { z = true; t &= 0x0FFFFFFF; }
    if (t > 1000 && t < 2000) { z = true; t -= 1000; }
    else if (t > 2000 && t < 3000) t -= 2000;  // M only
    else if (t > 3000) { z = true; t -= 3000; }
    auto point = [&]() { Vec3 v{dbl(), dbl(), z ? dbl() : 0.0}; return v; };
    auto line = [&]() { std::vector<Vec3> l; const std::uint32_t m = u32(); for (std::uint32_t k = 0; k < m; ++k) l.push_back(point()); return l; };
    switch (t) {
      case 1: f.type = 1; f.parts.push_back({point()}); break;
      case 2: f.type = 3; f.parts.push_back(line()); break;
      case 3: { f.type = 5; const std::uint32_t r = u32(); for (std::uint32_t k = 0; k < r; ++k) f.parts.push_back(line()); break; }
      case 4: case 5: case 6: case 7: { const std::uint32_t m = u32(); for (std::uint32_t k = 0; k < m; ++k) geom(f); break; }
      default: throw Error(ErrorCode::Unsupported, "WKB geometry type " + std::to_string(t));
    }
  }
};
}  // namespace

MapLayer read_geopackage(const std::filesystem::path& gpkg, const std::string& layer) {
#ifndef STRATONIX_HAVE_SQLITE
  (void)gpkg;
  (void)layer;
  throw Error(ErrorCode::Unsupported, "GeoPackage needs SQLite; this build has none (shapefiles and GeoJSON are read)");
#else
  MapLayer L;
  L.source = gpkg.string();
  L.name = layer;
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(gpkg.string().c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
    const std::string why = db ? sqlite3_errmsg(db) : "cannot open";
    sqlite3_close(db);
    throw Error(ErrorCode::Io, gpkg.string() + ": " + why);
  }
  std::string gcol;
  sqlite3_stmt* st = nullptr;
  if (sqlite3_prepare_v2(db, "SELECT column_name FROM gpkg_geometry_columns WHERE table_name = ?", -1, &st, nullptr) == SQLITE_OK) {
    sqlite3_bind_text(st, 1, layer.c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(st) == SQLITE_ROW) gcol = reinterpret_cast<const char*>(sqlite3_column_text(st, 0));
  }
  sqlite3_finalize(st);
  if (gcol.empty()) { sqlite3_close(db); throw Error(ErrorCode::Validation, gpkg.string() + ": no feature table '" + layer + "' in gpkg_geometry_columns"); }
  std::string q = "SELECT * FROM \"" + layer + "\"";
  if (sqlite3_prepare_v2(db, q.c_str(), -1, &st, nullptr) != SQLITE_OK) { sqlite3_close(db); throw Error(ErrorCode::Parse, gpkg.string() + ": cannot read " + layer); }
  int id = 0;
  while (sqlite3_step(st) == SQLITE_ROW) {
    MapFeature f;
    f.id = ++id;
    for (int c = 0; c < sqlite3_column_count(st); ++c) {
      const std::string name = sqlite3_column_name(st, c);
      if (name == gcol) {
        const auto* blob = static_cast<const unsigned char*>(sqlite3_column_blob(st, c));
        const std::size_t n = static_cast<std::size_t>(sqlite3_column_bytes(st, c));
        if (!blob || n < 8 || blob[0] != 'G' || blob[1] != 'P') continue;
        static const std::size_t env[8] = {0, 32, 48, 48, 64, 0, 0, 0};
        const std::size_t start = 8 + env[(blob[3] >> 1) & 7];
        if (blob[3] & 0x10) continue;  // empty geometry
        Wkb w{blob + start, n - start};
        try { w.geom(f); } catch (const Error& e) { L.report.push_back("feature " + std::to_string(id) + ": " + e.what()); }
      } else if (sqlite3_column_type(st, c) != SQLITE_NULL) {
        const auto* t = sqlite3_column_text(st, c);
        f.attributes[name] = t ? reinterpret_cast<const char*>(t) : "";
      }
    }
    L.features.push_back(std::move(f));
  }
  sqlite3_finalize(st);
  sqlite3_close(db);
  return L;
#endif
}

// ---- elevation models --------------------------------------------------------------------------
Grid2D read_esri_ascii_grid(const std::filesystem::path& file) {
  std::ifstream in(file);
  if (!in) throw Error(ErrorCode::Io, "cannot open " + file.string());
  std::size_t nc = 0, nr = 0;
  double xll = 0, yll = 0, cs = 0, nodata = std::numeric_limits<double>::quiet_NaN();
  bool centre = false;
  std::string key;
  std::streampos data = in.tellg();
  for (int k = 0; k < 6; ++k) {
    data = in.tellg();
    if (!(in >> key)) break;
    for (auto& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (key == "ncols") in >> nc;
    else if (key == "nrows") in >> nr;
    else if (key == "xllcorner" || key == "xllcenter") { in >> xll; centre = key == "xllcenter"; }
    else if (key == "yllcorner" || key == "yllcenter") in >> yll;
    else if (key == "cellsize") in >> cs;
    else if (key == "nodata_value") in >> nodata;
    else { in.seekg(data); break; }
  }
  if (!nc || !nr || cs <= 0) throw Error(ErrorCode::Parse, file.string() + ": an ESRI ASCII grid needs ncols, nrows and cellsize");
  const double off = centre ? 0.0 : 0.5 * cs;
  Grid2D g(GridSpec{xll + off, yll + off, cs, cs, nc, nr});
  for (std::size_t r = 0; r < nr; ++r)
    for (std::size_t c = 0; c < nc; ++c) {
      double v;
      if (!(in >> v)) throw Error(ErrorCode::Parse, file.string() + ": fewer values than ncols x nrows");
      g.at(c, nr - 1 - r) = (v == nodata) ? std::numeric_limits<double>::quiet_NaN() : v;  // the first row is the northern one
    }
  return g;
}

Grid2D read_geotiff_dem(const std::filesystem::path& file) {
  const auto b = read_bytes(file);
  if (b.size() < 8 || !((b[0] == 'I' && b[1] == 'I') || (b[0] == 'M' && b[1] == 'M'))) throw Error(ErrorCode::Parse, file.string() + ": not a TIFF");
  const bool le = b[0] == 'I';
  auto u16 = [&](std::size_t o) { return le ? unsigned(b[o] | (b[o + 1] << 8)) : unsigned((b[o] << 8) | b[o + 1]); };
  auto u32 = [&](std::size_t o) { return le ? le32(b.data() + o) : be32(b.data() + o); };
  auto f64 = [&](std::size_t o) { if (le) return led(b.data() + o); unsigned char t[8]; for (int k = 0; k < 8; ++k) t[k] = b[o + 7 - static_cast<std::size_t>(k)]; return led(t); };
  if (u16(2) != 42) throw Error(ErrorCode::Unsupported, file.string() + ": BigTIFF is not read");
  const std::size_t ifd = u32(4);
  const unsigned n = u16(ifd);
  struct Tag { unsigned type; std::uint32_t count; std::size_t off; };
  std::map<unsigned, Tag> tags;
  for (unsigned i = 0; i < n; ++i) {
    const std::size_t e = ifd + 2 + 12 * i;
    Tag t{u16(e + 2), u32(e + 4), 0};
    const std::size_t sz = (t.type == 3 ? 2 : t.type == 4 || t.type == 11 ? 4 : t.type == 12 ? 8 : 1) * t.count;
    t.off = sz <= 4 ? e + 8 : u32(e + 8);
    tags[u16(e)] = t;
  }
  auto val = [&](unsigned tag, std::size_t k = 0) -> double {
    auto it = tags.find(tag);
    if (it == tags.end()) return -1;
    const Tag& t = it->second;
    if (t.type == 3) return u16(t.off + 2 * k);
    if (t.type == 4) return u32(t.off + 4 * k);
    if (t.type == 12) return f64(t.off + 8 * k);
    return -1;
  };
  const std::size_t w = static_cast<std::size_t>(val(256)), h = static_cast<std::size_t>(val(257));
  const int bits = static_cast<int>(val(258)), comp = static_cast<int>(val(259)), fmt = tags.count(339) ? static_cast<int>(val(339)) : 1;
  if (comp != 1) throw Error(ErrorCode::Unsupported, file.string() + ": compressed GeoTIFF (compression " + std::to_string(comp) + ") is not read; write it uncompressed");
  if (tags.count(322)) throw Error(ErrorCode::Unsupported, file.string() + ": tiled GeoTIFF is not read; write it in strips");
  if (val(277) > 1) throw Error(ErrorCode::Unsupported, file.string() + ": more than one band");
  if (!tags.count(33550) || !tags.count(33922)) throw Error(ErrorCode::Validation, file.string() + ": no ModelPixelScale and ModelTiepoint: not georeferenced");
  const double sx = val(33550, 0), sy = val(33550, 1);
  const double I = val(33922, 0), J = val(33922, 1), X = val(33922, 3), Y = val(33922, 4);
  // pixel-is-area (default): tie point at the corner of the pixel; nodes at pixel centres
  bool point = false;
  if (tags.count(34735)) {
    const std::size_t nk = static_cast<std::size_t>(val(34735, 3));
    for (std::size_t k = 0; k < nk; ++k)
      if (val(34735, 4 + 4 * k) == 1025 && val(34735, 4 + 4 * k + 3) == 2) point = true;
  }
  const double half = point ? 0.0 : 0.5;
  const double x0 = X + (half - I) * sx, ytop = Y - (half - J) * sy;
  double nodata = std::numeric_limits<double>::quiet_NaN();
  if (tags.count(42113)) {
    const Tag& t = tags[42113];
    nodata = std::strtod(std::string(reinterpret_cast<const char*>(b.data() + t.off), t.count).c_str(), nullptr);
  }
  Grid2D g(GridSpec{x0, ytop - sy * static_cast<double>(h - 1), sx, sy, w, h});
  const std::size_t rps = tags.count(278) ? static_cast<std::size_t>(val(278)) : h, ns = tags[273].count;
  const int bps = bits / 8;
  for (std::size_t s = 0; s < ns; ++s) {
    std::size_t o = static_cast<std::size_t>(val(273, s));
    for (std::size_t r = s * rps; r < std::min(h, (s + 1) * rps); ++r)
      for (std::size_t c = 0; c < w; ++c, o += static_cast<std::size_t>(bps)) {
        if (o + static_cast<std::size_t>(bps) > b.size()) throw Error(ErrorCode::Parse, file.string() + ": strip beyond the end of the file");
        double v;
        if (fmt == 3 && bits == 32) { std::uint32_t x = u32(o); float f; std::memcpy(&f, &x, 4); v = f; }
        else if (fmt == 3 && bits == 64) v = f64(o);
        else if (bits == 16) v = fmt == 2 ? static_cast<double>(static_cast<std::int16_t>(u16(o))) : u16(o);
        else if (bits == 32) v = fmt == 2 ? static_cast<double>(static_cast<std::int32_t>(u32(o))) : u32(o);
        else if (bits == 8) v = fmt == 2 ? static_cast<double>(static_cast<std::int8_t>(b[o])) : b[o];
        else throw Error(ErrorCode::Unsupported, file.string() + ": " + std::to_string(bits) + "-bit samples");
        g.at(c, h - 1 - r) = (v == nodata) ? std::numeric_limits<double>::quiet_NaN() : v;  // the first row is the northern one
      }
  }
  return g;
}

// ---- map to model data -------------------------------------------------------------------------
namespace {
bool inside(const std::vector<std::vector<Vec3>>& rings, double x, double y) {
  bool in = false;  // even-odd over all rings (holes included)
  for (auto& r : rings)
    for (std::size_t i = 0, j = r.size() ? r.size() - 1 : 0; i < r.size(); j = i++)
      if (((r[i].y > y) != (r[j].y > y)) && (x < (r[j].x - r[i].x) * (y - r[i].y) / (r[j].y - r[i].y) + r[i].x)) in = !in;
  return in;
}
}  // namespace

ModelData map_to_model_data(const MapLayer& units, const std::string& field, const std::vector<std::string>& order, const MapLayer* faults,
                            const std::vector<Measurement>& ms, const Grid2D& dem) {
  ModelData d;
  const GridSpec& g = dem.spec;
  // unit at every elevation-model node
  std::vector<int> u(g.size(), -1);
  std::size_t unknown = 0;
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i)
      for (auto& f : units.features) {
        if (f.type != 5 || !inside(f.parts, g.x(i), g.y(j))) continue;
        auto it = f.attributes.find(field);
        const std::string name = it == f.attributes.end() ? "" : it->second;
        const auto k = std::find(order.begin(), order.end(), name);
        if (k == order.end()) { ++unknown; break; }
        u[g.index(i, j)] = static_cast<int>(k - order.begin());
        break;
      }
  if (unknown) d.report.push_back(std::to_string(unknown) + " elevation-model nodes lie in map units not in the stated unit order: left out");
  // contacts between neighbouring nodes of different units: a point on the base of the younger
  // unit, at the midpoint, on the elevation model
  std::size_t cross = 0;
  auto contact = [&](std::size_t a, std::size_t b) {
    const int ua = u[a], ub = u[b];
    if (ua < 0 || ub < 0 || ua == ub) return;
    const std::size_t young = ua < ub ? a : b, old = ua < ub ? b : a;
    const double x = 0.5 * (g.x(young % g.nx) + g.x(old % g.nx)), y = 0.5 * (g.y(young / g.nx) + g.y(old / g.nx));
    const double z = dem.sample(x, y);
    if (!std::isfinite(z)) return;
    if (std::abs(ua - ub) > 1) ++cross;
    d.contact_points.push_back({order[static_cast<std::size_t>(std::min(ua, ub))], {x, y, z}});
  };
  for (std::size_t j = 0; j < g.ny; ++j)
    for (std::size_t i = 0; i < g.nx; ++i) {
      if (i + 1 < g.nx) contact(g.index(i, j), g.index(i + 1, j));
      if (j + 1 < g.ny) contact(g.index(i, j), g.index(i, j + 1));
    }
  if (cross) d.report.push_back(std::to_string(cross) + " contact points between units not adjacent in the order (unconformities, faults or missing units)");
  for (auto& m : ms) {
    Measurement x = m;
    if (!std::isfinite(x.p.z) || x.p.z == 0) {
      const double z = dem.sample(x.p.x, x.p.y);
      if (std::isfinite(z)) x.p.z = z;
    }
    d.orientations.push_back(x);
  }
  if (faults)
    for (auto& f : faults->features)
      for (auto& part : f.parts) {
        std::vector<Vec3> tr;
        for (auto& p : part) { const double z = dem.sample(p.x, p.y); tr.push_back({p.x, p.y, std::isfinite(z) ? z : p.z}); }
        std::string nm = std::to_string(f.id);
        for (const char* k : {"name", "fname", "NAME", "FNAME", "id"}) { auto it = f.attributes.find(k); if (it != f.attributes.end() && !it->second.empty()) { nm = it->second; break; } }
        d.fault_traces.push_back({nm, tr});
      }
  d.report.push_back(std::to_string(d.contact_points.size()) + " contact points, " + std::to_string(d.orientations.size()) + " measurements, " +
                     std::to_string(d.fault_traces.size()) + " fault traces");
  return d;
}

namespace {
const desc::InputDesc& input_of(const desc::Description& d, const std::string& id) {
  for (auto& x : d.inputs)
    if (x.id == id) return x;
  throw Error(ErrorCode::Build, "input '" + id + "' is not listed under inputs");
}
MapLayer read_layer(const desc::Description& d, const desc::Description::MapSource& s) {
  const auto& in = input_of(d, s.input);
  MapLayer L = s.format == "shapefile" ? read_shapefile(in.resolved) : s.format == "geopackage" ? read_geopackage(in.resolved, s.layer) : read_geojson(in.resolved);
  // CRS: the input's stated CRS, else the file's; transformed to the model's when they differ (F46)
  CrsInfo from;
  bool known = false;
  if (in.has_crs) { from = in.crs; known = true; }
  else if (!L.crs_name.empty()) {
    const auto p = L.crs_name.rfind(':');
    if (p != std::string::npos) { from.epsg = std::atoi(L.crs_name.c_str() + p + 1); known = from.epsg > 0; }
  } else if (!L.crs_wkt.empty()) { from.wkt = L.crs_wkt; known = true; }
  if (!known) {
    L.report.push_back(s.input + ": no CRS stated by the input or the file: read as in the model's CRS");
    return L;
  }
  const bool same = (from.epsg > 0 && from.epsg == d.crs.epsg) || (!from.wkt.empty() && from.wkt == d.crs.wkt);
  if (same) { L.report.push_back(s.input + ": in the model's CRS (" + d.crs.label() + ")"); return L; }
  CrsTransform t(from, d.crs);
  std::vector<double> x, y, z;
  for (auto& f : L.features) for (auto& r : f.parts) for (auto& p : r) { x.push_back(p.x); y.push_back(p.y); z.push_back(p.z); }
  t.apply(x, y, z);
  std::size_t k = 0;
  for (auto& f : L.features) for (auto& r : f.parts) for (auto& p : r) { p.x = x[k]; p.y = y[k]; ++k; }
  L.report.push_back(s.input + ": transformed from " + from.label() + " by " + t.operation().name +
                     (t.operation().accuracy >= 0 ? " (stated accuracy " + std::to_string(t.operation().accuracy) + " m)" : " (no stated accuracy)"));
  return L;
}
}  // namespace

ModelData geological_map_data(const desc::Description& d) {
  const auto& G = d.geological_map;
  MapLayer units = read_layer(d, G.units);
  MapLayer faults;
  if (!G.faults.input.empty()) faults = read_layer(d, G.faults);
  std::vector<Measurement> ms;
  if (!G.measurements.input.empty()) {
    MapLayer m = read_layer(d, G.measurements);
    for (auto& f : m.features) {
      if (f.parts.empty() || f.parts[0].empty()) continue;
      if (!G.filter_field.empty()) {
        auto it = f.attributes.find(G.filter_field);
        if (it == f.attributes.end() || it->second != G.filter_value) continue;
      }
      auto a = f.attributes.find(G.dip_field), b = f.attributes.find(G.direction_field);
      if (a == f.attributes.end() || b == f.attributes.end()) continue;
      char* e1 = nullptr; char* e2 = nullptr;
      const double dip = std::strtod(a->second.c_str(), &e1), dir = std::strtod(b->second.c_str(), &e2);
      if (e1 == a->second.c_str() || e2 == b->second.c_str()) continue;
      Measurement x;
      x.p = f.parts[0][0];
      x.p.z = 0;
      x.dip = dip * kPi / 180.0;
      const double dd = G.convention == "strike_rhr" ? dir + 90.0 : dir;
      x.azimuth = std::fmod(std::fmod(dd, 360.0) + 360.0, 360.0) * kPi / 180.0;
      x.kind = G.filter_value.empty() ? "measurement" : G.filter_value;
      x.unit = "deg";
      x.source = G.source.empty() ? input_of(d, G.measurements.input).source : G.source;
      x.date = G.date;
      ms.push_back(x);
    }
  }
  const auto& em = input_of(d, G.dem_input);
  Grid2D dem = G.dem_format == "esri_ascii" ? read_esri_ascii_grid(em.resolved) : read_geotiff_dem(em.resolved);
  if (G.sample_spacing.present && G.sample_spacing.si() > dem.spec.dx) {
    // contacts sampled on a coarser lattice of the elevation model (bilinear)
    const double sp = G.sample_spacing.si();
    const double w = dem.spec.dx * static_cast<double>(dem.spec.nx - 1), h = dem.spec.dy * static_cast<double>(dem.spec.ny - 1);
    Grid2D c(GridSpec{dem.spec.x0, dem.spec.y0, sp, sp, static_cast<std::size_t>(w / sp) + 1, static_cast<std::size_t>(h / sp) + 1});
    for (std::size_t j = 0; j < c.spec.ny; ++j)
      for (std::size_t i = 0; i < c.spec.nx; ++i) c.at(i, j) = dem.sample(c.spec.x(i), c.spec.y(j));
    dem = c;
  }
  // map units with the same base form one model unit: their mutual contacts are internal. Units
  // without a base (the oldest, below the model's last horizon) form the basement group.
  std::vector<std::string> order;
  std::map<std::string, std::string> group_of;
  for (auto& [mu, base] : G.unit_order) {
    const std::string g = base.empty() ? std::string("(below the model)") : base;
    group_of[mu] = g;
    if (std::find(order.begin(), order.end(), g) == order.end()) order.push_back(g);
  }
  MapLayer grouped = units;
  for (auto& f : grouped.features) {
    auto it = f.attributes.find(G.units.field);
    if (it == f.attributes.end()) continue;
    auto g = group_of.find(it->second);
    it->second = g == group_of.end() ? "(not in the unit order) " + it->second : g->second;
  }
  ModelData md = map_to_model_data(grouped, G.units.field, order, G.faults.input.empty() ? nullptr : &faults, ms, dem);
  // each measurement belongs to the base horizon of the map unit it lies in
  std::size_t outside = 0;
  std::vector<Measurement> kept;
  for (auto& m : md.orientations) {
    for (auto& f : units.features) {
      if (f.type != 5 || !inside(f.parts, m.p.x, m.p.y)) continue;
      auto it = f.attributes.find(G.units.field);
      if (it == f.attributes.end()) break;
      for (auto& [mu, base] : G.unit_order)
        if (mu == it->second) m.horizon = base;  // empty for the basement group: left out
      break;
    }
    if (m.horizon.empty()) ++outside;
    else kept.push_back(m);
  }
  if (outside) md.report.push_back(std::to_string(outside) + " measurements lie outside the ordered map units: left out");
  // within the model's domain; repeated measurement locations kept once
  const double x0 = d.origin.si(0), y0 = d.origin.si(1), x1 = x0 + d.size.si(0), y1 = y0 + d.size.si(1);
  auto in_domain = [&](const Vec3& p) { return p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1; };
  std::vector<Measurement> uniq;
  std::size_t dup = 0, out_m = 0;
  for (auto& m : kept) {
    if (!in_domain(m.p)) { ++out_m; continue; }
    bool seen = false;
    for (auto& u : uniq) seen = seen || std::hypot(u.p.x - m.p.x, u.p.y - m.p.y) < 1.0;
    if (seen) { ++dup; continue; }
    uniq.push_back(m);
  }
  if (dup) md.report.push_back(std::to_string(dup) + " measurements repeat a location within 1 m: the first kept");
  if (out_m) md.report.push_back(std::to_string(out_m) + " measurements lie outside the model's domain: left out");
  md.orientations = uniq;
  const std::size_t nc = md.contact_points.size();
  md.contact_points.erase(std::remove_if(md.contact_points.begin(), md.contact_points.end(), [&](const auto& c) { return !in_domain(c.second); }),
                          md.contact_points.end());
  if (nc != md.contact_points.size()) md.report.push_back(std::to_string(nc - md.contact_points.size()) + " contact points outside the model's domain: left out");
  md.report.insert(md.report.begin(), units.report.begin(), units.report.end());
  md.report.push_back("elevation model " + std::to_string(dem.spec.nx) + " x " + std::to_string(dem.spec.ny) + " nodes at " + std::to_string(dem.spec.dx) + " m");
  return md;
}

}  // namespace stratonix::io
