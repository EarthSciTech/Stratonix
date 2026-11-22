#include "stratonix/core/well.hpp"

#include <cmath>

#include "stratonix/core/error.hpp"

namespace stratonix {

WellPath::WellPath(std::string name, double x, double y, double datum_z, double total_depth)
    : name_(std::move(name)), x_(x), y_(y), datum_z_(datum_z) {
  pts_ = {{x, y, datum_z}, {x, y, -total_depth}};
  md_ = {0.0, total_depth + datum_z};
}

WellPath::WellPath(std::string name, double x, double y, double datum_z, std::vector<WellStation> st)
    : name_(std::move(name)), x_(x), y_(y), datum_z_(datum_z), stations_(std::move(st)) {
  if (stations_.size() < 2) throw Error(ErrorCode::InvalidArgument, "a deviated well needs at least two stations");
  Vec3 p{x, y, datum_z};
  pts_.push_back(p);
  md_.push_back(stations_[0].md);
  for (std::size_t i = 1; i < stations_.size(); ++i) {
    const auto& a = stations_[i - 1];
    const auto& b = stations_[i];
    double dmd = b.md - a.md;
    if (dmd <= 0) throw Error(ErrorCode::InvalidArgument, "well stations must increase in measured depth");
    // Minimum curvature.
    double cosb = std::cos(b.inc - a.inc) - std::sin(a.inc) * std::sin(b.inc) * (1 - std::cos(b.azi - a.azi));
    double beta = std::acos(std::clamp(cosb, -1.0, 1.0));
    double rf = beta > 1e-12 ? 2.0 / beta * std::tan(beta / 2) : 1.0;
    double dn = dmd / 2 * (std::sin(a.inc) * std::cos(a.azi) + std::sin(b.inc) * std::cos(b.azi)) * rf;
    double de = dmd / 2 * (std::sin(a.inc) * std::sin(a.azi) + std::sin(b.inc) * std::sin(b.azi)) * rf;
    double dv = dmd / 2 * (std::cos(a.inc) + std::cos(b.inc)) * rf;
    p = {p.x + de, p.y + dn, p.z - dv};
    pts_.push_back(p);
    md_.push_back(b.md);
  }
}

WellPath WellPath::from_positions(std::string name, double x, double y, double datum_z, const std::vector<double>& md,
                                  const std::vector<Vec3>& net, std::string method) {
  if (md.size() < 2 || md.size() != net.size()) throw Error(ErrorCode::InvalidArgument, "a deviated well needs at least two positions");
  WellPath w(std::move(name), x, y, datum_z, md.back() - md.front());
  w.pts_.clear();
  w.md_ = md;
  for (auto& p : net) w.pts_.push_back({x + p.y, y + p.x, datum_z - p.z});
  w.stations_ = {{md.front(), 0, 0}, {md.back(), 0, 0}};  // marks the path as deviated
  w.method_ = std::move(method);
  return w;
}

Vec3 WellPath::at_md(double m) const {
  if (pts_.empty()) return {x_, y_, datum_z_ - m};
  if (md_.size() < 2 || m <= md_.front()) return {pts_.front().x, pts_.front().y, pts_.front().z - (m - (md_.empty() ? 0 : md_.front()))};
  for (std::size_t i = 1; i < md_.size(); ++i)
    if (m <= md_[i]) {
      const double t = (m - md_[i - 1]) / (md_[i] - md_[i - 1]);
      return pts_[i - 1] + (pts_[i] - pts_[i - 1]) * t;
    }
  return {pts_.back().x, pts_.back().y, pts_.back().z - (m - md_.back())};
}

Vec3 WellPath::at_elevation(double z) const {
  if (pts_.empty()) return {x_, y_, z};
  if (z >= pts_.front().z) return {pts_.front().x, pts_.front().y, z};
  for (std::size_t i = 1; i < pts_.size(); ++i) {
    const auto& a = pts_[i - 1];
    const auto& b = pts_[i];
    if (z <= a.z && z >= b.z) {
      double t = a.z == b.z ? 0 : (a.z - z) / (a.z - b.z);
      return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), z};
    }
  }
  // Below total depth: extend the last segment's direction vertically.
  return {pts_.back().x, pts_.back().y, z};
}

double WellPath::md_at_elevation(double z) const {
  if (pts_.size() < 2) return datum_z_ - z;
  for (std::size_t i = 1; i < pts_.size(); ++i) {
    const auto& a = pts_[i - 1];
    const auto& b = pts_[i];
    if (z <= a.z && z >= b.z) {
      double t = a.z == b.z ? 0 : (a.z - z) / (a.z - b.z);
      return md_[i - 1] + t * (md_[i] - md_[i - 1]);
    }
  }
  return md_.back() + (pts_.back().z - z);
}

}  // namespace stratonix
