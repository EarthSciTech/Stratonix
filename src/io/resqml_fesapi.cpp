// RESQML 2.0.1 export of a corner-point grid through FESAPI (F2I-Consulting, Apache 2.0), built
// only with STRATONIX_WITH_FESAPI (STX-TN-12 section 6). One EPC package and its HDF5 file:
// a local depth CRS, an IjkGridExplicitRepresentation with its coordinate lines and split
// coordinate lines (faults), the active cells, and ZONE, BLOCK and the grid properties.
#include <cmath>
#include <map>
#include <tuple>
#include <type_traits>

#include <fesapi/common/DataObjectRepository.h>
#include <fesapi/common/EpcDocument.h>
#include <fesapi/eml2/AbstractHdfProxy.h>
#include <fesapi/resqml2/ContinuousProperty.h>
#include <fesapi/resqml2/DiscreteProperty.h>
#include <fesapi/resqml2/IjkGridExplicitRepresentation.h>
#if __has_include(<fesapi/eml2/AbstractLocal3dCrs.h>)
#include <fesapi/eml2/AbstractLocal3dCrs.h>
#endif
#if __has_include(<fesapi/resqml2/AbstractLocal3dCrs.h>)
#include <fesapi/resqml2/AbstractLocal3dCrs.h>
#endif

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/io/export.hpp"

namespace stratonix::io {

namespace fs = std::filesystem;

namespace {
// Element types of FESAPI's array arguments differ between releases (unsigned int or uint64_t);
// they are read from the member function's own signature.
template <class T> struct Args;
template <class C, class R, class... A> struct Args<R (C::*)(A...)> { using type = std::tuple<A...>; };
using GeomArgs = Args<decltype(&RESQML2_NS::IjkGridExplicitRepresentation::setGeometryAsCoordinateLineNodes)>::type;
template <std::size_t N> using ElemOf = std::remove_cv_t<std::remove_pointer_t<std::tuple_element_t<N, GeomArgs>>>;
using PillarIndex = ElemOf<6>;   // pillarOfCoordinateLine
using CumCount = ElemOf<7>;      // splitCoordinateLineColumnCumulativeCount
using ColumnIndex = ElemOf<8>;   // splitCoordinateLineColumns
}  // namespace

namespace {
auto* repr_active(COMMON_NS::DataObjectRepository& repo, RESQML2_NS::IjkGridExplicitRepresentation* ijk) {
  return repo.createDiscreteProperty(ijk, "", "ACTIVE", 1, gsoap_eml2_3::eml23__IndexableElement::cells,
                                     gsoap_resqml2_0_1::resqml20__ResqmlPropertyKind::active);
}
}  // namespace

std::vector<ExportedFile> write_resqml_grid(ExportContext& c, const fs::path& folder) {
  const auto& g = *c.cells.grid;
  const int nx = g.nx, ny = g.ny, nz = g.nz;
  const std::size_t np = static_cast<std::size_t>(nx + 1) * (ny + 1);
  // Nodes per coordinate line: the top of every layer and the base of the last; a layer whose base
  // is not the next layer's top (a gap between zones) takes the next layer's top (reported).
  auto node_depth = [&](int i, int j, int k, int corner) {  // depth of interface k at corner (0..3) of column (i, j)
    if (k < nz) return -g.zcorn[g.zcorn_index(i, j, k, corner)];
    return -g.zcorn[g.zcorn_index(i, j, nz - 1, corner + 4)];
  };
  double gap = 0;
  for (int k = 0; k + 1 < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        for (int q = 0; q < 4; ++q)
          gap = std::max(gap, std::fabs(g.zcorn[g.zcorn_index(i, j, k, q + 4)] - g.zcorn[g.zcorn_index(i, j, k + 1, q)]));
  // Coordinate lines: per pillar, the columns around it grouped by identical node depths; the
  // first group uses the pillar's own line, the others are split coordinate lines.
  struct Line { std::size_t pillar; std::vector<std::size_t> columns; std::vector<double> depth; };
  std::vector<Line> lines(np);
  std::vector<Line> splits;
  for (int pj = 0; pj <= ny; ++pj)
    for (int pi = 0; pi <= nx; ++pi) {
      const std::size_t p = g.pillar(pi, pj);
      std::vector<Line> groups;
      for (int q = 0; q < 4; ++q) {
        int ci = pi - 1 + (q & 1), cj = pj - 1 + (q >> 1);
        if (ci < 0 || cj < 0 || ci >= nx || cj >= ny) continue;
        const int corner = (q & 1 ? 0 : 1) + (q >> 1 ? 0 : 2);  // the column's corner at this pillar
        std::vector<double> dz(static_cast<std::size_t>(nz) + 1);
        for (int k = 0; k <= nz; ++k) dz[static_cast<std::size_t>(k)] = node_depth(ci, cj, k, corner);
        bool placed = false;
        for (auto& gr : groups) {
          bool same = true;
          for (std::size_t k = 0; k < dz.size() && same; ++k) same = std::fabs(gr.depth[k] - dz[k]) < 1e-6;
          if (same) { gr.columns.push_back(static_cast<std::size_t>(ci + nx * cj)); placed = true; break; }
        }
        if (!placed) groups.push_back({p, {static_cast<std::size_t>(ci + nx * cj)}, dz});
      }
      lines[p] = groups.front();
      for (std::size_t s = 1; s < groups.size(); ++s) splits.push_back(groups[s]);
    }
  const std::size_t nlines = np + splits.size();
  std::vector<double> pts(3 * nlines * (static_cast<std::size_t>(nz) + 1));
  auto put_node = [&](std::size_t line, int k, std::size_t p, double depth) {
    const double* cd = &g.coord[6 * p];
    const double dt = -cd[2], db = -cd[5];
    const double t = std::fabs(db - dt) > 1e-12 ? (depth - dt) / (db - dt) : 0.0;
    double* o = &pts[3 * (static_cast<std::size_t>(k) * nlines + line)];
    o[0] = cd[0] + t * (cd[3] - cd[0]);
    o[1] = cd[1] + t * (cd[4] - cd[1]);
    o[2] = depth;
  };
  for (int k = 0; k <= nz; ++k) {
    for (std::size_t p = 0; p < np; ++p) put_node(p, k, p, lines[p].depth[static_cast<std::size_t>(k)]);
    for (std::size_t s = 0; s < splits.size(); ++s) put_node(np + s, k, splits[s].pillar, splits[s].depth[static_cast<std::size_t>(k)]);
  }
  std::vector<PillarIndex> pillar_of(splits.size());
  std::vector<CumCount> cum(splits.size());
  std::vector<ColumnIndex> cols;
  for (std::size_t s = 0; s < splits.size(); ++s) {
    pillar_of[s] = static_cast<PillarIndex>(splits[s].pillar);
    for (auto col : splits[s].columns) cols.push_back(static_cast<ColumnIndex>(col));
    cum[s] = static_cast<CumCount>(cols.size());
  }

  fs::create_directories(folder / "grid");
  const std::string base = petrel_name(c.info().name) + "_grid";
  const fs::path epc = folder / "grid" / (base + ".epc");
  COMMON_NS::DataObjectRepository repo;
  EML2_NS::AbstractHdfProxy* hdf = repo.createHdfProxy("", "Hdf Proxy", (folder / "grid").string(), base + ".h5",
                                                       COMMON_NS::DataObjectRepository::openingMode::OVERWRITE);
  repo.setDefaultHdfProxy(hdf);
  auto* crs = repo.createLocalDepth3dCrs("", "Stratonix local depth CRS", 0.0, 0.0, 0.0, 0.0, gsoap_resqml2_0_1::eml20__LengthUom::m,
                                         c.info().crs.name.empty() ? std::string("local") : c.info().crs.name, gsoap_resqml2_0_1::eml20__LengthUom::m,
                                         "depth positive down from the model datum", false);
  repo.setDefaultCrs(crs);
  auto* ijk = repo.createIjkGridExplicitRepresentation("", c.info().name + " corner-point grid", static_cast<unsigned int>(nx), static_cast<unsigned int>(ny),
                                                       static_cast<unsigned int>(nz));
  ijk->setGeometryAsCoordinateLineNodes(gsoap_resqml2_0_1::resqml20__PillarShape::straight, gsoap_resqml2_0_1::resqml20__KDirection::down, false,
                                        pts.data(), hdf, splits.size(), pillar_of.data(), cum.data(), cols.data(), nullptr, crs);
  // active cells: a discrete property of the RESQML kind "active" (read as such by resqpy)
  {
    std::vector<std::int8_t> act(g.ncell());
    for (std::size_t n = 0; n < act.size(); ++n) act[n] = static_cast<std::int8_t>(g.actnum[n] != 0);
    auto* ap = repr_active(repo, ijk);
    ap->pushBackInt8Hdf5Array3dOfValues(act.data(), static_cast<std::uint64_t>(nx), static_cast<std::uint64_t>(ny), static_cast<std::uint64_t>(nz), hdf, -1);
  }
  std::size_t nprops = 0;
  if (c.cells.props)
    for (auto& p : c.cells.props->props) {
      if (p.integer) {
        std::vector<std::int64_t> v(p.v.size());
        for (std::size_t n = 0; n < v.size(); ++n) v[n] = std::isnan(p.v[n]) ? -1 : std::llround(p.v[n]);
        auto* dp = repo.createDiscreteProperty(ijk, "", p.name, 1, gsoap_eml2_3::eml23__IndexableElement::cells,
                                               gsoap_resqml2_0_1::resqml20__ResqmlPropertyKind::discrete);
        dp->pushBackInt64Hdf5Array3dOfValues(v.data(), static_cast<std::uint64_t>(nx), static_cast<std::uint64_t>(ny), static_cast<std::uint64_t>(nz), hdf, -1);
      } else {
        auto* cp = repo.createContinuousProperty(ijk, "", p.name, 1, gsoap_eml2_3::eml23__IndexableElement::cells,
                                                 gsoap_resqml2_0_1::resqml20__ResqmlUom::Euc, gsoap_resqml2_0_1::resqml20__ResqmlPropertyKind::continuous);
        std::vector<double> v = p.v;  // NaN is the missing value
        cp->pushBackDoubleHdf5Array3dOfValues(v.data(), static_cast<std::uint64_t>(nx), static_cast<std::uint64_t>(ny), static_cast<std::uint64_t>(nz), hdf);
      }
      ++nprops;
    }
  COMMON_NS::EpcDocument doc(epc.string());
  doc.serializeFrom(repo);
  hdf->close();
  return {{"grid/" + base + ".epc", "RESQML 2.0.1 EPC package (FESAPI)", "RESQML 2.0.1 (*.epc)",
           "IJK grid " + std::to_string(nx) + " x " + std::to_string(ny) + " x " + std::to_string(nz) + ", " + std::to_string(splits.size()) +
               " split coordinate lines (faults), " + std::to_string(nprops) + " properties; layer gaps closed by up to " + fmt_fixed(gap, 3) + " m"},
          {"grid/" + base + ".h5", "HDF5 (RESQML arrays)", "-", "arrays of the EPC package above"}};
}

}  // namespace stratonix::io
