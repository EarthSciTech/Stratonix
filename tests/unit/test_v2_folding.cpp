// Stage 2B (Prompts 20 and 21): trishear, fault-bend folds by fault-parallel flow, detachment
// folds, parallel and similar mechanics, en echelon arrays (STX-TN-06).
#include "v2_fixture.hpp"

using namespace v2fx;

namespace {

struct P2 { double u = 0, z = 0; };  // dip-plane coordinates: u east (dip direction), z elevation
P2 operator+(P2 a, P2 b) { return {a.u + b.u, a.z + b.z}; }
P2 operator-(P2 a, P2 b) { return {a.u - b.u, a.z - b.z}; }
P2 operator*(P2 a, double s) { return {a.u * s, a.z * s}; }
double dot2(P2 a, P2 b) { return a.u * b.u + a.z * b.z; }

/// Independent trishear (Zehnder and Allmendinger 2000, symmetric linear field) for a fault through
/// (0, -1300) dipping east; v0 = +1 moves the hanging wall up dip (reverse).
struct Trishear {
  double dip, phi, ps, slip, v0;
  P2 tip1;  // present tip
  P2 ex() const { return {-std::cos(dip), std::sin(dip)}; }   // up dip, towards the region ahead of the tip
  P2 ny() const { return {std::sin(dip), std::cos(dip)}; }    // normal into the hanging wall
  P2 tip_at(double s) const { return tip1 - ex() * (ps * (slip - s)); }
  P2 vel(P2 p, double s) const {
    P2 d = p - tip_at(s);
    double x = dot2(d, ex()), y = dot2(d, ny()), tp = std::tan(phi);
    double vx, vy;
    if (x <= 0) { vx = y > 0 ? v0 : 0; vy = 0; }
    else if (y >= x * tp) { vx = v0; vy = 0; }
    else if (y <= -x * tp) { vx = 0; vy = 0; }
    else {
      double r = y / (x * tp);
      vx = 0.5 * v0 * (r + 1);
      vy = 0.25 * v0 * tp * (r * r - 1);
    }
    return ex() * vx + ny() * vy;
  }
  P2 forward(P2 p) const {
    const int n = 4000;
    double h = slip / n;
    for (int i = 0; i < n; ++i) {
      double s = i * h;
      P2 k1 = vel(p, s), k2 = vel(p + k1 * (h / 2), s + h / 2), k3 = vel(p + k2 * (h / 2), s + h / 2), k4 = vel(p + k3 * h, s + h);
      p = p + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);
    }
    return p;
  }
};

std::string trishear_fault(const std::string& dip, const std::string& tip_depth, const std::string& throw_keys,
                           const std::string& fold, const std::string& upper = "") {
  return "  - name: F1\n    age: " + q("1 ordinal") + "\n    geometry: {type: planar, reference_point: " + q("[5000, 5000] m") +
         ", reference_depth: " + q("1300 m") + ", strike: " + q("0 deg") + ", dip: " + q(dip) +
         "}\n    extent: {upper_tip: " + (upper.empty() ? "{rule: depth, depth: " + q(tip_depth) + "}" : upper) + ", lower_tip: " + kFree + ", length: " +
         q("100000 m") + "}\n    throw: {" + throw_keys + ", profile_strike: " + flat_top("0.01") +
         "}\n    slip_mode: dip_slip\n    hanging_wall: {type: rigid}\n    fold: " + fold + "\n";
}
const std::string kTri = "{type: trishear, apical_angle: " + q("60 deg") + ", propagation_slip_ratio: " + q("1.5 ratio") + "}";

double scalar_at(const Model& m, int h, double x, double z) {
  EvalResult r;
  std::vector<Vec3> p{{x, 5000, z}};
  m.evaluate(p, kEvalScalars, r);
  return r.s(static_cast<std::size_t>(h), 0);
}

// Reverse trishear: slip 400 m on a 30-degree fault, present tip at 1200 m.
Trishear reverse_case() {
  Trishear t{deg2rad(30), deg2rad(30), 1.5, 400, 1.0, {0, 0}};
  t.tip1 = {-100 / std::tan(t.dip), -1200};
  return t;
}

}  // namespace

TEST_CASE("P20 trishear: every horizon against an independent integration") {
  auto m = build_text(model_yaml(trishear_fault("30 deg", "1200 m", "max: " + q("200 m") + ", sense: reverse", kTri)));
  Trishear t = reverse_case();
  const double depths[] = {1000, 1200, 1400, 1600};
  double worst = 0;
  int n = 0;
  for (int h = 0; h < 4; ++h)
    for (double u = -2500; u <= 2500; u += 50) {
      P2 p = t.forward({u, -depths[h]});
      // skip material within 30 m of the fault plane behind the tip (cut-offs, overlap)
      P2 d = p - t.tip1;
      if (dot2(d, t.ex()) < 0 && std::fabs(dot2(d, t.ny())) < 30) continue;
      worst = std::max(worst, std::fabs(scalar_at(*m, h, 5000 + p.u, p.z)));
      ++n;
    }
  std::cout << "  " << n << " points, worst horizon misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(n > 300);
  CHECK(worst <= 0.5);
}

TEST_CASE("P20 trishear round trip and area") {
  desc::Description d;
  auto m = build_text(model_yaml(trishear_fault("30 deg", "1200 m", "max: " + q("200 m") + ", sense: reverse", kTri)), &d);
  Trishear t = reverse_case();
  double worst = 0, worst_det = 0;
  int n = 0;
  for (double u = -1500; u <= 1500; u += 100)
    for (double z = -700; z >= -1700; z -= 100) {
      P2 p0{u, z}, p = t.forward(p0);
      std::vector<Vec3> in{{5000 + p.u, 5000, p.z}}, out(1);
      m->restore(in, out);
      worst = std::max(worst, std::hypot(out[0].x - 5000 - p0.u, out[0].z - p0.z));
      // Jacobian of the restoration in the dip plane (central differences, 0.5 m)
      const double e = 0.5;
      std::vector<Vec3> q4{{5000 + p.u + e, 5000, p.z}, {5000 + p.u - e, 5000, p.z}, {5000 + p.u, 5000, p.z + e}, {5000 + p.u, 5000, p.z - e}}, r4(4);
      m->restore(q4, r4);
      double a = (r4[0].x - r4[1].x) / (2 * e), b = (r4[2].x - r4[3].x) / (2 * e), c = (r4[0].z - r4[1].z) / (2 * e), dd = (r4[2].z - r4[3].z) / (2 * e);
      worst_det = std::max(worst_det, std::fabs(a * dd - b * c - 1));
      ++n;
    }
  std::cout << "  " << n << " points: worst round trip " << fmt_fixed(worst, 4) << " m, worst |det - 1| " << fmt_fixed(worst_det, 5) << "\n";
  CHECK(worst <= 0.05);
  CHECK(worst_det <= 0.01);
  d.wells.clear();
  auto rep = checks_for(*m, d);
  auto* c12 = find_check(rep, "C12", "F1");
  REQUIRE(c12 != nullptr);
  CHECK(c12->passed);
}

TEST_CASE("P20 forced fold above a blind normal fault") {
  auto m = build_text(model_yaml(trishear_fault("60 deg", "1300 m", "max: " + q("100 m") + ", sense: normal",
                                                "{type: trishear, apical_angle: " + q("60 deg") + ", propagation_slip_ratio: " + q("1 ratio") + "}")));
  // Far from the fold: footwall unmoved, hanging wall lowered by the throw.
  CHECK_NEAR(depth_faulted(*m, 0, 2000, 5000), 1000, 0.01);
  CHECK_NEAR(depth_faulted(*m, 0, 8000, 5000), 1100, 0.05);
  // H1 lies above the tip: folded, not cut, so it is continuous and never rises eastwards.
  double prev = 0;
  bool ok = true;
  for (double x = 2000; x <= 8000; x += 25) {
    double z = depth_faulted(*m, 0, x, 5000);
    ok = ok && !std::isnan(z) && z >= prev - 1e-6;
    prev = z;
  }
  CHECK(ok);
}

namespace {
/// Independent fault-parallel flow on a polyline (u, depth) profile; slip > 0 moves the hanging wall
/// towards smaller u (up dip, reverse).
struct FPF {
  std::vector<P2> pts;  // (u, depth)
  P2 e(std::size_t k) const { P2 d = pts[k + 1] - pts[k]; double l = std::hypot(d.u, d.z); return d * (1 / l); }
  P2 n(std::size_t k) const { P2 t = e(k); return {t.z, -t.u}; }  // towards the hanging wall (up)
  P2 miter(std::size_t k) const {
    if (k == 0) return n(0);
    if (k + 1 == pts.size()) return n(k - 1);
    P2 a = n(k - 1), b = n(k);
    return (a + b) * (1 / (1 + dot2(a, b)));
  }
  P2 V(std::size_t k, double d) const { return pts[k] + miter(k) * d; }
  // (s, d) of a hanging-wall point
  void coords(P2 p, double& s, double& d) const {
    std::size_t ns = pts.size() - 1;
    for (std::size_t k = 0; k < ns; ++k) {
      d = dot2(p - pts[k], n(k));
      double t = dot2(p - V(k, d), e(k));
      double len = dot2(V(k + 1, d) - V(k, d), e(k));
      if ((k == 0 || t >= 0) && (k + 1 == ns || t <= len)) {
        s = t;
        for (std::size_t j = 0; j < k; ++j) s += dot2(V(j + 1, d) - V(j, d), e(j));
        return;
      }
    }
  }
  P2 at(double s, double d) const {
    std::size_t ns = pts.size() - 1;
    double acc = 0;
    for (std::size_t k = 0; k < ns; ++k) {
      double len = dot2(V(k + 1, d) - V(k, d), e(k));
      if (s < acc + len || k + 1 == ns) return V(k, d) + e(k) * (s - acc);
      acc += len;
    }
    return {};
  }
  P2 forward(P2 p, double slip) const { double s, d; coords(p, s, d); return at(s - slip, d); }
  // hanging wall: above the fault (the footwall does not move)
  bool above(P2 p) const { double s, d; coords(p, s, d); return d > 0; }
};
// Thrust: upper flat at 1100 m, 30-degree ramp to a lower flat at 2000 m.
const std::vector<P2> kThrust = {{-6000, 1100}, {2000, 1100}, {2000 + 900 / std::tan(deg2rad(30)), 2000}, {12000, 2000}};
std::string fpf_fault() {
  return "  - name: F1\n    age: " + q("1 ordinal") + "\n    geometry: {type: ramp_flat, reference_point: " + q("[5000, 5000] m") +
         ", strike: " + q("0 deg") + ", profile: " +
         q("[[-6000, 1100], [2000, 1100], [" + fmt_fixed(2000 + 900 / std::tan(deg2rad(30)), 6) + ", 2000], [12000, 2000]] m") +
         "}\n    extent: {upper_tip: " + kFree + ", lower_tip: " + kFree + ", length: " + q("100000 m") +
         "}\n    throw: {slip: " + q("1000 m") + ", sense: reverse, profile_strike: " + flat_top("0.01") +
         "}\n    slip_mode: dip_slip\n    hanging_wall: {type: fault_parallel_flow}\n";
}
}  // namespace

TEST_CASE("P21 fault-bend fold by fault-parallel flow") {
  auto m = build_text(model_yaml(fpf_fault()));
  FPF f{kThrust};
  // Hanging-wall horizons against the independent flow (the reference point is at u = 0, x = 5000).
  double worst = 0;
  int n = 0;
  for (int h = 2; h < 4; ++h)
    for (double u = 2600; u <= 9000; u += 50) {
      double z0 = 1400 + 200 * (h - 2);
      if (!f.above({u, z0})) continue;
      P2 p = f.forward({u, z0}, 1000);
      if (p.u > 9000) continue;
      worst = std::max(worst, std::fabs(scalar_at(*m, h, 5000 + p.u, -p.z)));
      ++n;
    }
  std::cout << "  " << n << " points, worst misfit " << fmt_fixed(worst, 4) << " m\n";
  CHECK(worst <= 0.5);
  // H3 above the ramp: parallel to it (30 degrees) at its original 600 m above the lower flat.
  double z1 = depth_in_block(*m, 2, 1, 5000 + 3100, 5000), z2 = depth_in_block(*m, 2, 1, 5000 + 3400, 5000);
  CHECK_NEAR(rad2deg(std::atan((z2 - z1) / 300)), 30, 0.5);
  double nrm = (3100 - 2000) * std::sin(deg2rad(30)) - (z1 - 1100) * std::cos(deg2rad(30));  // distance above the ramp line
  CHECK_NEAR(nrm, 600, 0.5);
  // Far behind the ramp the hanging wall only moved along the flat.
  CHECK_NEAR(depth_in_block(*m, 2, 1, 5000 + 7000, 5000), 1400, 0.01);
}

TEST_CASE("P21 fault-parallel flow round trip") {
  auto m = build_text(model_yaml(fpf_fault()));
  FPF f{kThrust};
  double worst = 0;
  for (double u = 2700; u <= 8000; u += 150)
    for (double z = 1300; z <= 1900; z += 100) {
      if (!f.above({u, z})) continue;
      P2 p = f.forward({u, z}, 1000);
      std::vector<Vec3> in{{5000 + p.u, 5000, -p.z}}, out(1);
      m->restore(in, out);
      worst = std::max(worst, std::hypot(out[0].x - 5000 - u, out[0].z + z));
    }
  std::cout << "  worst round trip " << fmt_fixed(worst, 5) << " m\n";
  CHECK(worst <= 0.01);
}

namespace {
std::string prim(const std::string& body) { return "structure:\n  regional: {type: flat}\n  primitives:\n    - " + body + "\n"; }
}  // namespace

TEST_CASE("P21 detachment fold: uplifted area equals shortening times depth to the detachment") {
  auto m = build_text(model_yaml("", prim("{type: detachment_fold, name: D, crest: " + q("[5000, 5000] m") + ", axis_azimuth: " +
                                          q("0 deg") + ", half_width: " + q("2000 m") + ", shortening: " + q("300 m") +
                                          ", detachment: H4}")));
  const double depths[] = {1000, 1200, 1400, 1600};
  for (int h = 0; h < 4; ++h) {
    double area = 0;
    for (double x = 2000; x <= 8000; x += 5) area += (depths[h] - depth_faulted(*m, h, x, 5000)) * 5;
    CHECK_NEAR(area, 300 * (1600 - depths[h]), 0.01 * 300 * 600 + 1);
  }
}

TEST_CASE("P21 parallel and similar mechanics") {
  std::string a = "{type: anticline, name: A, crest: " + q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") + ", half_width_left: " +
                  q("1500 m") + ", half_width_right: " + q("1500 m") + ", relief: {H1: " + q("300 m") + "}, mechanics: ";
  {
    auto m = build_text(model_yaml("", prim(a + "parallel}")));
    // normal distance from H2 to H1 equals the regional thickness (200 m)
    std::vector<P2> h1;
    for (double x = 1000; x <= 9000; x += 2) h1.push_back({x, -depth_faulted(*m, 0, x, 5000)});
    double worst = 0;
    for (double x = 3000; x <= 7000; x += 100) {
      P2 p{x, -depth_faulted(*m, 1, x, 5000)};
      double best = 1e9;
      for (auto& q1 : h1) best = std::min(best, std::hypot(q1.u - p.u, q1.z - p.z));
      worst = std::max(worst, std::fabs(best - 200));
    }
    std::cout << "  parallel: worst normal-thickness error " << fmt_fixed(worst, 3) << " m\n";
    CHECK(worst <= 2.0);
  }
  {
    auto m = build_text(model_yaml("", prim(a + "similar}")));
    const double depths[] = {1000, 1200, 1400, 1600};
    for (int h = 0; h < 4; ++h) CHECK_NEAR(depth_faulted(*m, h, 5000, 5000), depths[h] - 300, 0.01);
  }
}

TEST_CASE("P21 en echelon fold array") {
  auto m = build_text(model_yaml("", prim("{type: anticline, name: A, crest: " + q("[2000, 5000] m") + ", axis_azimuth: " + q("0 deg") +
                                          ", half_length: " + q("2500 m") + ", half_width_left: " + q("1000 m") + ", half_width_right: " +
                                          q("1000 m") + ", relief: {H1: " + q("100 m") + "}, array: {count: 3, step_along: " + q("500 m") +
                                          ", step_across: " + q("3000 m") + "}}")));
  int crests = 0;
  double a = depth_faulted(*m, 0, 0, 5250), b = depth_faulted(*m, 0, 25, 5250);
  for (double x = 50; x <= 10000; x += 25) {
    double c = depth_faulted(*m, 0, x, 5250);
    if (b < a - 1e-9 && b < c - 1e-9) ++crests;
    a = b;
    b = c;
  }
  CHECK(crests == 3);
}

TEST_CASE("P21 refusals (F17, F19)") {
  auto r1 = validate(model_yaml(trishear_fault("30 deg", "1200 m", "max: " + q("200 m"), kTri, kFree)));
  CHECK(!r1.ok());
  CHECK(r1.text().find("F17") != std::string::npos);
  auto r2 = validate(model_yaml("", prim("{type: detachment_fold, name: D, crest: " + q("[5000, 5000] m") + ", axis_azimuth: " + q("0 deg") +
                                         ", half_width: " + q("2000 m") + ", shortening: " + q("300 m") + ", detachment: H1}")));
  CHECK(!r2.ok());
  CHECK(r2.text().find("F19") != std::string::npos);
}

STX_TEST_MAIN()
