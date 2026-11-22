// Uncertainty (stage 3C, STX-TN-13): seeded ensembles of a description, key-number statistics,
// sensitivity and probability maps.
#include "stratonix/package/uncertainty.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <thread>

#include "stratonix/checks/checks.hpp"
#include "stratonix/checks/manifest.hpp"
#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"
#include "stratonix/engines/posterior.hpp"
#include "stratonix/engines/wells.hpp"
#include "stratonix/repr/topology.hpp"
#include "stratonix/figures/figures.hpp"
#include "stratonix/io/export.hpp"
#include "stratonix/io/io.hpp"
#include "stratonix/package/pipeline.hpp"
#include "stratonix/repr/repr.hpp"

namespace stratonix {
// Defined in pipeline.cpp.
fig::FigureContext figure_context(const desc::Description& d, const Model& m, const std::vector<std::pair<std::string, std::string>>& keys);
}  // namespace stratonix

namespace stratonix::unc {

namespace fs = std::filesystem;

// ---- generator ---------------------------------------------------------------------------------
namespace {
std::uint64_t splitmix64(std::uint64_t& x) {
  std::uint64_t z = (x += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
inline std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
}  // namespace

std::uint64_t member_seed(std::uint64_t seed, std::uint64_t member) {
  std::uint64_t x = seed ^ (member * 0xD1B54A32D192ED03ull);
  splitmix64(x);
  return splitmix64(x);
}
Rng::Rng(std::uint64_t seed) {
  std::uint64_t x = seed;
  for (auto& s : s_) s = splitmix64(x);
}
std::uint64_t Rng::next() {
  const std::uint64_t r = rotl(s_[1] * 5, 7) * 9, t = s_[1] << 17;
  s_[2] ^= s_[0];
  s_[3] ^= s_[1];
  s_[1] ^= s_[2];
  s_[0] ^= s_[3];
  s_[2] ^= t;
  s_[3] = rotl(s_[3], 45);
  return r;
}
double Rng::uniform() { return (static_cast<double>(next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }

// ---- distributions -------------------------------------------------------------------------------
double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
double normal_quantile(double p) {
  if (!(p > 0 && p < 1)) {
    if (p == 0) return -1e300;
    if (p == 1) return 1e300;
    return kNaN;
  }
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  const double plow = 0.02425;
  double x;
  if (p < plow) {
    double q = std::sqrt(-2 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= 1 - plow) {
    double q = p - 0.5, r = q * q;
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    double q = std::sqrt(-2 * std::log(1 - p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  // one Halley step on the CDF
  const double e = normal_cdf(x) - p, u = e * std::sqrt(2 * kPi) * std::exp(x * x / 2);
  return x - u / (1 + x * u / 2);
}

double Distribution::quantile(double u) const {
  if (type == "uniform") return a + u * (b - a);
  if (type == "normal") {
    const double fa = lo > -1e299 ? normal_cdf((lo - a) / b) : 0.0, fb = hi < 1e299 ? normal_cdf((hi - a) / b) : 1.0;
    return a + b * normal_quantile(fa + u * (fb - fa));
  }
  if (type == "triangular") {  // a = min, b = mode, c = max
    const double fc = (b - a) / (c - a);
    return u < fc ? a + std::sqrt(u * (c - a) * (b - a)) : c - std::sqrt((1 - u) * (c - a) * (c - b));
  }
  if (type == "discrete" || type == "posterior") {
    double cum = 0;
    for (std::size_t k = 0; k < values.size(); ++k) {
      cum += weights[k];
      if (u <= cum + 1e-15) return values[k];
    }
    return values.back();
  }
  throw std::invalid_argument("unknown distribution '" + type + "'");
}
double Distribution::mean() const {
  if (type == "uniform") return 0.5 * (a + b);
  if (type == "normal") {
    if (lo < -1e299 && hi > 1e299) return a;
    const double al = (lo - a) / b, be = (hi - a) / b, Z = normal_cdf(be) - normal_cdf(al);
    auto phi = [](double x) { return std::exp(-0.5 * x * x) / std::sqrt(2 * kPi); };
    return a + b * ((lo > -1e299 ? phi(al) : 0) - (hi < 1e299 ? phi(be) : 0)) / Z;
  }
  if (type == "triangular") return (a + b + c) / 3;
  double m = 0;
  for (std::size_t k = 0; k < values.size(); ++k) m += weights[k] * values[k];
  return m;
}
double Distribution::variance() const {
  if (type == "uniform") return (b - a) * (b - a) / 12;
  if (type == "normal") {
    if (lo < -1e299 && hi > 1e299) return b * b;
    const double al = (lo - a) / b, be = (hi - a) / b, Z = normal_cdf(be) - normal_cdf(al);
    auto phi = [](double x) { return std::exp(-0.5 * x * x) / std::sqrt(2 * kPi); };
    const double pa = lo > -1e299 ? phi(al) : 0, pb = hi < 1e299 ? phi(be) : 0;
    const double ta = lo > -1e299 ? al * pa : 0, tb = hi < 1e299 ? be * pb : 0;
    return b * b * (1 + (ta - tb) / Z - ((pa - pb) / Z) * ((pa - pb) / Z));
  }
  if (type == "triangular") return (a * a + b * b + c * c - a * b - a * c - b * c) / 18;
  const double m = mean();
  double v = 0;
  for (std::size_t k = 0; k < values.size(); ++k) v += weights[k] * (values[k] - m) * (values[k] - m);
  return v;
}

// ---- specification -------------------------------------------------------------------------------
EnsembleSpec ensemble_spec(const desc::Description& d) {
  const auto& u = d.uncertainty;
  if (!u.present) throw Error(ErrorCode::InvalidArgument, "the description has no uncertainty block (STX-TN-13 section 6)", "uncertainty");
  EnsembleSpec s;
  s.members = static_cast<int>(std::llround(u.members.si()));
  s.seed = static_cast<std::uint64_t>(std::llround(u.seed.si()));
  s.one_at_a_time = u.one_at_a_time;
  s.key_numbers = u.key_numbers;
  for (auto& q : u.parameters) {
    const desc::Param* base = d.param_by_name(q.path);
    if (!base) throw Error(ErrorCode::Validation, "uncertainty: no parameter " + q.path, "uncertainty");
    UncertainParameter p;
    p.path = q.path;
    p.unit = base->unit;
    p.base = base->values.at(0);
    const double to = base->to_si;  // distribution parameters into the base parameter's unit
    auto v = [&](const desc::Param& x, std::size_t i = 0) { return x.si(i) / to; };
    p.dist.type = q.distribution;
    if (q.distribution == "uniform") { p.dist.a = v(q.min); p.dist.b = v(q.max); }
    if (q.distribution == "normal") {
      p.dist.a = v(q.mean);
      p.dist.b = v(q.sd);
      if (q.min.present) p.dist.lo = v(q.min);
      if (q.max.present) p.dist.hi = v(q.max);
    }
    if (q.distribution == "triangular") { p.dist.a = v(q.min); p.dist.b = v(q.mode); p.dist.c = v(q.max); }
    if (q.distribution == "discrete")
      for (std::size_t i = 0; i < q.values.values.size(); ++i) {
        p.dist.values.push_back(v(q.values, i));
        p.dist.weights.push_back(q.weights.si(i));
      }
    for (const desc::Param* x : {&q.mean, &q.sd, &q.min, &q.max, &q.mode, &q.values, &q.weights})
      if (x->present) {
        std::string ev = desc::evidence_name(x->evidence);
        if (p.evidence.find(ev) == std::string::npos) p.evidence += (p.evidence.empty() ? "" : " / ") + ev;
        if (!x->source.empty() && p.source.find(x->source) == std::string::npos) p.source += (p.source.empty() ? "" : "; ") + x->source;
      }
    s.parameters.push_back(std::move(p));
  }
  auto index_of = [&](const std::string& path) {
    for (std::size_t k = 0; k < s.parameters.size(); ++k)
      if (s.parameters[k].path == path) return static_cast<int>(k);
    return -1;
  };
  for (auto& c : u.correlations) s.correlations.push_back({{index_of(c.a), index_of(c.b)}, c.rho.si()});
  for (auto& sc : u.scenarios) {
    Scenario x;
    x.name = sc.name;
    for (auto& o : sc.options) x.options.push_back({o.name, o.weight.si(), o.set, o.remove});
    s.scenarios.push_back(std::move(x));
  }
  for (auto& m : u.maps)
    s.maps.push_back({m.kind, m.unit, m.horizon, m.depth.present ? m.depth.si() : 0.0, m.spacing.present ? m.spacing.si() : 100.0});
  if (!s.correlations.empty()) {
    std::vector<double> corr(s.parameters.size() * s.parameters.size(), 0.0);
    const int n = static_cast<int>(s.parameters.size());
    for (int i = 0; i < n; ++i) corr[static_cast<std::size_t>(i * n + i)] = 1;
    for (auto& [ij, rho] : s.correlations) {
      corr[static_cast<std::size_t>(ij.first * n + ij.second)] = rho;
      corr[static_cast<std::size_t>(ij.second * n + ij.first)] = rho;
    }
    cholesky(corr, n);  // refuses a matrix that is not positive definite
  }
  return s;
}

std::vector<double> cholesky(const std::vector<double>& A, int n) {
  std::vector<double> L(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = A[static_cast<std::size_t>(i * n + j)];
      for (int k = 0; k < j; ++k) s -= L[static_cast<std::size_t>(i * n + k)] * L[static_cast<std::size_t>(j * n + k)];
      if (i == j) {
        if (!(s > 1e-12)) throw Error(ErrorCode::Validation, "F33: the correlation matrix is not positive definite", "uncertainty.correlations");
        L[static_cast<std::size_t>(i * n + i)] = std::sqrt(s);
      } else {
        L[static_cast<std::size_t>(i * n + j)] = s / L[static_cast<std::size_t>(j * n + j)];
      }
    }
  return L;
}

std::vector<double> correlated_uniforms(const EnsembleSpec& s, Rng& rng) {
  const int n = static_cast<int>(s.parameters.size());
  std::vector<double> u(static_cast<std::size_t>(n));
  for (auto& x : u) x = rng.uniform();
  if (s.correlations.empty()) return u;
  std::vector<double> corr(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) corr[static_cast<std::size_t>(i * n + i)] = 1;
  for (auto& [ij, rho] : s.correlations) {
    corr[static_cast<std::size_t>(ij.first * n + ij.second)] = rho;
    corr[static_cast<std::size_t>(ij.second * n + ij.first)] = rho;
  }
  auto L = cholesky(corr, n);
  std::vector<double> z(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = normal_quantile(u[static_cast<std::size_t>(i)]);
  for (int i = 0; i < n; ++i) {
    double y = 0;
    for (int k = 0; k <= i; ++k) y += L[static_cast<std::size_t>(i * n + k)] * z[static_cast<std::size_t>(k)];
    u[static_cast<std::size_t>(i)] = std::clamp(normal_cdf(y), 1e-16, 1 - 1e-16);
  }
  return u;
}

std::vector<MemberSample> sample_members(const EnsembleSpec& s) {
  std::vector<MemberSample> out;
  MemberSample base;
  base.index = 0;
  for (auto& p : s.parameters) base.values.push_back(p.base);
  base.options.assign(s.scenarios.size(), 0);
  out.push_back(base);
  for (int k = 1; k <= s.members; ++k) {
    MemberSample m;
    m.index = k;
    m.seed = member_seed(s.seed, static_cast<std::uint64_t>(k));
    Rng rng(m.seed);
    auto u = correlated_uniforms(s, rng);
    for (std::size_t i = 0; i < s.parameters.size(); ++i) m.values.push_back(s.parameters[i].dist.quantile(u[i]));
    for (auto& sc : s.scenarios) {
      const double v = rng.uniform();
      double cum = 0;
      int pick = static_cast<int>(sc.options.size()) - 1;
      for (std::size_t o = 0; o < sc.options.size(); ++o) {
        cum += sc.options[o].weight;
        if (v <= cum) { pick = static_cast<int>(o); break; }
      }
      m.options.push_back(pick);
    }
    out.push_back(std::move(m));
  }
  return out;
}

// ---- description documents -------------------------------------------------------------------
namespace {
std::vector<std::string> split_path(const std::string& path) {
  std::vector<std::string> out;
  std::size_t a = 0;
  while (a <= path.size()) {
    auto b = path.find('.', a);
    out.push_back(path.substr(a, b == std::string::npos ? std::string::npos : b - a));
    if (b == std::string::npos) break;
    a = b + 1;
  }
  return out;
}
/// Node at a path: map keys, and list items by their name. With `parent`, also the list and
/// index of a list item (for removal).
desc::Node* find_node(desc::Node& root, const std::string& path, desc::Node** list = nullptr, std::size_t* index = nullptr) {
  desc::Node* n = &root;
  for (auto& seg : split_path(path)) {
    if (n->is_seq()) {
      desc::Node* hit = nullptr;
      for (std::size_t k = 0; k < n->seq.size(); ++k) {
        const desc::Node* nm = n->seq[k].get("name");
        if (nm && nm->str == seg) {
          if (list) *list = n;
          if (index) *index = k;
          hit = &n->seq[k];
          break;
        }
      }
      n = hit;
    } else {
      if (list) *list = nullptr;
      n = n->get(seg);
    }
    if (!n) return nullptr;
  }
  return n;
}
std::string scalar_text(const desc::Node& n) { return n.is_number() ? fmt_shortest(n.number) : n.str; }
}  // namespace

bool get_parameter(const desc::Node& doc, const std::string& path, double* value, std::string* unit) {
  auto* n = find_node(const_cast<desc::Node&>(doc), path);
  if (!n) return false;
  std::string text;
  if (n->is_map()) {
    const desc::Node *v = n->get("value"), *u = n->get("unit");
    if (!v) return false;
    if (value) *value = v->is_number() ? v->number : std::stod(v->str);
    if (unit) *unit = u ? u->str : "";
    return true;
  }
  text = scalar_text(*n);
  auto bar = text.find('|');
  std::string vu = text.substr(0, bar);
  // "value unit": the last token is the unit
  while (!vu.empty() && vu.back() == ' ') vu.pop_back();
  auto sp = vu.rfind(' ');
  if (sp == std::string::npos) return false;
  try {
    if (value) *value = std::stod(vu.substr(0, sp));
  } catch (...) {
    return false;
  }
  if (unit) *unit = vu.substr(sp + 1);
  return true;
}

void set_parameter(desc::Node& doc, const std::string& path, double value, const std::string& unit, const std::string& evidence,
                   const std::string& source) {
  auto* n = find_node(doc, path);
  if (!n) throw Error(ErrorCode::InvalidArgument, "no parameter at '" + path + "'", path);
  std::string src = source;
  for (auto& ch : src)
    if (ch == '|') ch = '/';
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.10g", value);
  *n = desc::Node::make_string(std::string(buf) + " " + unit + " | " + evidence + " | " + src);
}

desc::Node member_document(const desc::Node& base, const EnsembleSpec& s, const MemberSample& m) {
  desc::Node doc = base;
  // the members carry no uncertainty block of their own; posterior members are not refitted
  const bool drop_cond = s.posterior && m.index != 0;
  if (doc.is_map())
    doc.map.erase(std::remove_if(doc.map.begin(), doc.map.end(),
                                 [&](auto& kv) { return kv.first == "uncertainty" || (drop_cond && kv.first == "conditioning"); }),
                  doc.map.end());
  if (m.index != 0)  // members (> 0) and one-at-a-time runs (< 0); the base case keeps its values
    for (std::size_t i = 0; i < s.parameters.size(); ++i) {
      const auto& p = s.parameters[i];
      const std::string d = p.dist.type == "uniform"      ? "uniform(" + fmt_shortest(p.dist.a) + ", " + fmt_shortest(p.dist.b) + ")"
                            : p.dist.type == "normal"     ? "normal(" + fmt_shortest(p.dist.a) + ", " + fmt_shortest(p.dist.b) + ")"
                            : p.dist.type == "triangular" ? "triangular(" + fmt_shortest(p.dist.a) + ", " + fmt_shortest(p.dist.b) + ", " + fmt_shortest(p.dist.c) + ")"
                                                          : "discrete";
      if (p.dist.type == "posterior") {
        set_parameter(doc, p.path, m.values[i], p.unit, "calibrated",
                      p.source + (m.index > 0 ? "; member " + std::to_string(m.index) : std::string("; one-at-a-time run")));
        continue;
      }
      set_parameter(doc, p.path, m.values[i], p.unit, "assumed",
                    "sampled: " + d + " " + p.unit + " (" + p.evidence + (p.source.empty() ? "" : ", " + p.source) + "; base " + fmt_shortest(p.base) +
                        " " + p.unit + ")" + (m.index > 0 ? "; member " + std::to_string(m.index) + ", seed " + std::to_string(m.seed) : std::string("; one-at-a-time run")));
    }
  for (std::size_t k = 0; k < s.scenarios.size(); ++k) {
    const auto& o = s.scenarios[k].options[static_cast<std::size_t>(m.options[k])];
    for (auto& [path, node] : o.set) {
      auto segs = split_path(path);
      std::string parent_path, leaf = segs.back();
      for (std::size_t q = 0; q + 1 < segs.size(); ++q) parent_path += (q ? "." : "") + segs[q];
      desc::Node* parent = parent_path.empty() ? &doc : find_node(doc, parent_path);
      if (!parent) throw Error(ErrorCode::Validation, "scenario " + s.scenarios[k].name + "/" + o.name + ": no node at '" + parent_path + "'", "uncertainty");
      parent->set(leaf, node);
    }
    for (auto& path : o.remove) {
      desc::Node* list = nullptr;
      std::size_t idx = 0;
      if (!find_node(doc, path, &list, &idx) || !list)
        throw Error(ErrorCode::Validation, "scenario " + s.scenarios[k].name + "/" + o.name + ": no list item at '" + path + "'", "uncertainty");
      list->seq.erase(list->seq.begin() + static_cast<std::ptrdiff_t>(idx));
    }
  }
  return doc;
}

// ---- statistics ------------------------------------------------------------------------------------
Stats statistics(std::vector<double> v) {
  Stats s;
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
  s.n = v.size();
  if (v.empty()) {
    s.mean = s.sd = s.p10 = s.p50 = s.p90 = s.min = s.max = kNaN;
    return s;
  }
  std::sort(v.begin(), v.end());
  double sum = 0;
  for (double x : v) sum += x;
  s.mean = sum / static_cast<double>(v.size());
  double ss = 0;
  for (double x : v) ss += (x - s.mean) * (x - s.mean);
  s.sd = v.size() > 1 ? std::sqrt(ss / static_cast<double>(v.size() - 1)) : 0.0;
  auto pct = [&](double p) {
    const double h = (static_cast<double>(v.size()) - 1) * p;
    const std::size_t lo = static_cast<std::size_t>(std::floor(h));
    const std::size_t hi = std::min(v.size() - 1, lo + 1);
    return v[lo] + (h - static_cast<double>(lo)) * (v[hi] - v[lo]);
  };
  s.p10 = pct(0.1);
  s.p50 = pct(0.5);
  s.p90 = pct(0.9);
  s.min = v.front();
  s.max = v.back();
  return s;
}

double first_order_index(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = std::min(x.size(), y.size());
  if (n < 4) return kNaN;
  std::vector<std::size_t> ord(n);
  for (std::size_t i = 0; i < n; ++i) ord[i] = i;
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
  double my = 0;
  for (std::size_t i = 0; i < n; ++i) my += y[i];
  my /= static_cast<double>(n);
  double vy = 0;
  for (std::size_t i = 0; i < n; ++i) vy += (y[i] - my) * (y[i] - my);
  vy /= static_cast<double>(n);
  if (!(vy > 0)) return 0.0;
  // a factor with few distinct values (a scenario's option): the exact means per value
  {
    std::map<double, std::pair<double, std::size_t>> groups;
    for (std::size_t i = 0; i < n; ++i) {
      auto& g = groups[x[i]];
      g.first += y[i];
      ++g.second;
    }
    if (groups.size() <= 10) {
      double vg = 0;
      for (auto& [v, g] : groups) {
        const double m = g.first / static_cast<double>(g.second);
        vg += static_cast<double>(g.second) / static_cast<double>(n) * (m - my) * (m - my);
      }
      return vg / vy;
    }
  }
  const std::size_t bins = std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(n)))));
  double vb = 0;
  for (std::size_t b = 0; b < bins; ++b) {
    const std::size_t lo = b * n / bins, hi = (b + 1) * n / bins;
    if (hi <= lo) continue;
    double m = 0;
    for (std::size_t i = lo; i < hi; ++i) m += y[ord[i]];
    m /= static_cast<double>(hi - lo);
    vb += static_cast<double>(hi - lo) / static_cast<double>(n) * (m - my) * (m - my);
  }
  return vb / vy;
}

std::size_t EnsembleResult::ok() const {
  std::size_t n = 0;
  for (auto& m : members) n += m.sample.index > 0 && m.status == "ok";
  return n;
}

bool parse_key_value(const std::string& text, double* value, std::string* unit) {
  std::size_t i = 0;
  while (i < text.size() && text[i] == ' ') ++i;
  std::size_t j = i;
  if (j < text.size() && (text[j] == '-' || text[j] == '+')) ++j;
  bool digit = false;
  while (j < text.size() && (std::isdigit(static_cast<unsigned char>(text[j])) || text[j] == '.')) { digit = digit || text[j] != '.'; ++j; }
  if (!digit) return false;
  if (value) *value = std::stod(text.substr(i, j - i));
  std::size_t k = j;
  while (k < text.size() && text[k] == ' ') ++k;
  std::size_t e = k;
  while (e < text.size() && text[e] != ' ' && text[e] != '(' && text[e] != ',') ++e;
  if (unit) *unit = text.substr(k, e - k);
  return true;
}

// ---- the runner -------------------------------------------------------------------------------------
namespace {

struct MemberRun {
  MemberResult r;
  std::unique_ptr<ContactGraph> graph;  ///< stage 4A: for the topology comparison
  std::vector<std::vector<std::uint8_t>> map_hits;  ///< per probability map, per node (1 hit)
};

bool reported_key(const std::string& k) { return k != "Evidence basis" && k != "Horizons / faults" && k != "Checks"; }

MemberRun run_member(const desc::Node& base, const fs::path& base_dir, const EnsembleSpec& s, const MemberSample& m,
                     const std::vector<GridSpec>& map_grids) {
  MemberRun out;
  out.r.sample = m;
  out.map_hits.resize(s.maps.size());
  desc::Description d;
  try {
    auto doc = member_document(base, s, m);
    auto rep = desc::validate_description_text(desc::to_json(doc), base_dir, &d);
    if (!rep.ok()) {
      out.r.status = "refused";
      std::string why;
      for (auto& is : rep.issues)
        if (is.error) { why = is.rule + " " + is.path + ": " + is.message; break; }
      out.r.reason = why;
      return out;
    }
  } catch (const std::exception& e) {
    out.r.status = "refused";
    out.r.reason = e.what();
    return out;
  }
  ModelPtr model;
  try {
    BuildOptions o;
    o.figures = false;
    model = build_only(d, o).model;
  } catch (const std::exception& e) {
    out.r.status = "build_failed";
    out.r.reason = e.what();
    return out;
  }
  try {
    CheckContext ctx(*model, d, make_well_paths(d), d.grid_spacing.si());
    auto rep = run_checks(ctx);
    out.graph = std::make_unique<ContactGraph>(ctx.contact_graph());
    for (auto& c : rep.results)
      if (!c.passed && !c.warn_only && !c.info) out.r.failed_checks.push_back(c.id + " " + c.subject);
    for (auto& [k, v] : key_numbers(d, *model, rep, ctx)) {
      double x;
      std::string u;
      if (reported_key(k) && parse_key_value(v, &x, &u)) out.r.keys.push_back({k, x});
    }
    out.r.status = rep.passed() ? "ok" : "checks_failed";
    if (!rep.passed()) out.r.reason = std::to_string(out.r.failed_checks.size()) + " checks failed";
  } catch (const std::exception& e) {
    out.r.status = "build_failed";
    out.r.reason = std::string("checks: ") + e.what();
    return out;
  }
  // probability maps
  for (std::size_t k = 0; k < s.maps.size(); ++k) {
    const auto& spec = s.maps[k];
    const auto& g = map_grids[k];
    auto& hits = out.map_hits[k];
    hits.assign(g.size(), 0);
    if (spec.kind == "unit") {
      int unit = -1;
      for (std::size_t u = 0; u < model->info().units.size(); ++u)
        if (model->info().units[u].name == spec.unit) unit = static_cast<int>(u);
      std::vector<Vec3> pts(g.size());
      for (std::size_t j = 0; j < g.ny; ++j)
        for (std::size_t i = 0; i < g.nx; ++i) pts[g.index(i, j)] = {g.x(i), g.y(j), -spec.depth};
      EvalResult er;
      model->evaluate(pts, kEvalUnit, er);
      for (std::size_t n = 0; n < pts.size(); ++n) hits[n] = er.unit[n] == unit;
    } else {
      int h = -1;
      for (std::size_t q = 0; q < model->info().horizons.size(); ++q)
        if (model->info().horizons[q].name == spec.horizon) h = static_cast<int>(q);
      auto fsurf = repr::horizon_grid_faulted(*model, h, g);
      const bool seal = d.closure_faults_seal.present && d.closure_faults_seal.flag;
      auto cl = repr::closure_analysis(fsurf.z, seal);
      for (std::size_t n = 0; n < hits.size(); ++n) hits[n] = cl.mask.v[n] == 1.0;
    }
  }
  return out;
}

std::string csv_text(const std::string& s) {
  if (s.find_first_of(",\"\n") == std::string::npos) return s;
  std::string o = "\"";
  for (char c : s) o += c == '"' ? std::string("\"\"") : std::string(1, c);
  return o + "\"";
}

}  // namespace

EnsembleResult run_ensemble(const fs::path& description, const fs::path& out_dir, const EnsembleOptions& o) {
  register_all_plugins();
  auto d = desc::load_description(description);
  EnsembleResult res;
  const bool bayes = d.conditioning.present && d.conditioning.method == "bayesian";
  res.spec = (d.uncertainty.present || !bayes) ? ensemble_spec(d) : EnsembleSpec{};
  auto& s = res.spec;
  if (o.members > 0) s.members = o.members;
  if (o.seed >= 0) s.seed = static_cast<std::uint64_t>(o.seed);
  // stage 4A (STX-TN-14 section 4): a Bayesian conditioning gives the members of its posterior,
  // one per retained sample (every thin-th step of every walker), or an even subset of them
  PosteriorResult post;
  std::vector<std::size_t> post_rows;
  if (bayes) {
    BuildInputs in;
    io::load_build_inputs(d, in);
    BuildLog lg;
    post = sample_posterior(d, in, posterior_options(d), lg);
    if (!post.ok) throw Error(ErrorCode::Build, "posterior sampling refused: " + post.reason, "conditioning");
    std::vector<std::size_t> rows;
    for (std::size_t w = 0; w < static_cast<std::size_t>(post.walkers); ++w)
      for (std::size_t st = static_cast<std::size_t>(post.thin) - 1; st < static_cast<std::size_t>(post.steps); st += static_cast<std::size_t>(post.thin))
        rows.push_back(w * static_cast<std::size_t>(post.steps) + st);
    const std::size_t n = (o.members > 0 && static_cast<std::size_t>(o.members) < rows.size()) ? static_cast<std::size_t>(o.members) : rows.size();
    for (std::size_t k = 0; k < n; ++k) post_rows.push_back(rows[k * rows.size() / n]);
    s.members = static_cast<int>(n);
    if (!d.uncertainty.present) s.seed = post.seed;
    // the conditioned parameters leave the prior sampling (their priors are inside the posterior)
    std::vector<int> keep;
    std::vector<UncertainParameter> kept;
    for (std::size_t i = 0; i < s.parameters.size(); ++i) {
      bool cond = std::find(post.paths.begin(), post.paths.end(), s.parameters[i].path) != post.paths.end();
      keep.push_back(cond ? -1 : static_cast<int>(kept.size()));
      if (!cond) kept.push_back(s.parameters[i]);
    }
    std::vector<std::pair<std::pair<int, int>, double>> corr;
    for (auto& [ij, rho] : s.correlations)
      if (keep[static_cast<std::size_t>(ij.first)] >= 0 && keep[static_cast<std::size_t>(ij.second)] >= 0)
        corr.push_back({{keep[static_cast<std::size_t>(ij.first)], keep[static_cast<std::size_t>(ij.second)]}, rho});
    s.parameters = kept;
    s.correlations = corr;
  }
  // every scenario option must give a valid description on its own
  for (std::size_t k = 0; k < s.scenarios.size(); ++k)
    for (std::size_t q = 0; q < s.scenarios[k].options.size(); ++q) {
      MemberSample t;
      for (auto& p : s.parameters) t.values.push_back(p.base);
      t.options.assign(s.scenarios.size(), 0);
      t.options[k] = static_cast<int>(q);
      desc::Description tmp;
      auto rep = desc::validate_description_text(desc::to_json(member_document(d.resolved, s, t)), description.parent_path(), &tmp);
      if (!rep.ok())
        throw Error(ErrorCode::Validation, "uncertainty: scenario " + s.scenarios[k].name + ", option " + s.scenarios[k].options[q].name +
                                               " gives an invalid description:\n" + rep.text(), "uncertainty.scenarios");
    }
  // the base case and the domain: the least-squares fit (the posterior is not sampled again)
  desc::Node base_doc = d.resolved;
  desc::Description dbase = d;
  if (bayes) {
    if (desc::Node* cn = find_node(base_doc, "conditioning"); cn && cn->is_map())
      cn->map.erase(std::remove_if(cn->map.begin(), cn->map.end(), [](auto& kv) { return kv.first == "method" || kv.first == "sampler"; }), cn->map.end());
    dbase.conditioning.method = "least_squares";
    dbase.conditioning.walkers = dbase.conditioning.steps = dbase.conditioning.burn_in = dbase.conditioning.thin = dbase.conditioning.seed = desc::Param{};
  }
  auto samples = sample_members(s);
  if (bayes) {
    s.posterior = true;
    for (std::size_t k = 0; k < post.paths.size(); ++k) {
      const desc::Param* bp = d.param_by_name(post.paths[k]);
      UncertainParameter up;
      up.path = post.paths[k];
      up.unit = bp->unit;
      up.base = bp->values.at(0);
      up.evidence = "calibrated";
      up.source = "posterior sample of conditioning.method: bayesian (" + std::to_string(post.walkers) + " walkers, seed " + std::to_string(post.seed) + ")";
      up.dist.type = "posterior";
      for (std::size_t r : post_rows) up.dist.values.push_back(post.chain[r][k] / bp->to_si);
      std::sort(up.dist.values.begin(), up.dist.values.end());
      up.dist.weights.assign(up.dist.values.size(), 1.0 / static_cast<double>(up.dist.values.size()));
      samples[0].values.push_back(up.base);
      for (std::size_t m = 1; m < samples.size(); ++m) samples[m].values.push_back(post.chain[post_rows[m - 1]][k] / bp->to_si);
      s.parameters.push_back(std::move(up));
    }
  }
  // one-at-a-time runs: each parameter at its P10 and P90, the others at their base values
  std::vector<MemberSample> oat;
  if (s.one_at_a_time)
    for (std::size_t i = 0; i < s.parameters.size(); ++i)
      for (double q : {0.1, 0.9}) {
        MemberSample m = samples[0];
        m.index = -1;  // not an ensemble member: no seed, never in the statistics
        m.values[i] = s.parameters[i].dist.quantile(q);
        oat.push_back(m);
      }
  const Domain dom = build_only(dbase, BuildOptions{}).model->info().domain;
  std::vector<GridSpec> map_grids;
  for (auto& mp : s.maps) {
    GridSpec g;
    g.dx = g.dy = mp.spacing;
    g.nx = static_cast<std::size_t>(std::max(1.0, std::floor(dom.sx / mp.spacing)));
    g.ny = static_cast<std::size_t>(std::max(1.0, std::floor(dom.sy / mp.spacing)));
    g.x0 = dom.x0 + 0.5 * (dom.sx - static_cast<double>(g.nx - 1) * mp.spacing);
    g.y0 = dom.y0 + 0.5 * (dom.sy - static_cast<double>(g.ny - 1) * mp.spacing);
    map_grids.push_back(g);
  }
  // run every member (and the one-at-a-time runs) in parallel; results kept in order
  std::vector<MemberSample> all = samples;
  all.insert(all.end(), oat.begin(), oat.end());
  std::vector<MemberRun> runs(all.size());
  const unsigned saved_threads = thread_count();
  const unsigned nt = std::max(1u, std::min<unsigned>(o.threads ? o.threads : saved_threads, static_cast<unsigned>(all.size())));
  set_thread_count(1);  // members in parallel, each single-threaded inside
  const LogLevel saved_level = log_level();
  set_log_level(LogLevel::Error);
  std::atomic<std::size_t> next{0};
  auto worker = [&] {
    for (;;) {
      std::size_t k = next.fetch_add(1);
      if (k >= all.size()) return;
      runs[k] = run_member(base_doc, description.parent_path(), s, all[k], map_grids);
    }
  };
  {
    std::vector<std::thread> ts;
    for (unsigned t = 1; t < nt; ++t) ts.emplace_back(worker);
    worker();
    for (auto& t : ts) t.join();
  }
  set_log_level(saved_level);
  set_thread_count(saved_threads);
  for (std::size_t k = 0; k < samples.size(); ++k) res.members.push_back(runs[k].r);
  // stage 4A (STX-TN-14 section 6): each member's contact graph against the base case's
  if (runs[0].graph)
    for (std::size_t k = 1; k < samples.size(); ++k) {
      if (!runs[k].graph) continue;
      std::vector<std::string> diff;
      if (same_topology(*runs[0].graph, *runs[k].graph, &diff)) {
        res.members[k].topology = "same";
      } else {
        ++res.topology_differs;
        std::string t;
        for (std::size_t q = 0; q < diff.size() && q < 4; ++q) t += (q ? "; " : "") + diff[q];
        res.members[k].topology = t + (diff.size() > 4 ? "; and " + std::to_string(diff.size() - 4) + " more" : "");
      }
    }

  // key names in the order of first appearance (base case first)
  std::vector<std::string> keys;
  for (auto& m : res.members)
    for (auto& [k, v] : m.keys)
      if (std::find(keys.begin(), keys.end(), k) == keys.end()) keys.push_back(k);
  if (!s.key_numbers.empty()) {
    std::vector<std::string> sel;
    for (auto& k : keys)
      if (std::find(s.key_numbers.begin(), s.key_numbers.end(), k) != s.key_numbers.end()) sel.push_back(k);
    keys = sel;
  }
  auto value_of = [](const MemberResult& m, const std::string& k) {
    for (auto& [n, v] : m.keys)
      if (n == k) return v;
    return kNaN;
  };
  std::size_t failed = 0;
  for (auto& m : res.members) failed += m.sample.index > 0 && m.status != "ok";
  for (auto& k : keys) {
    std::vector<double> v;
    for (auto& m : res.members)
      if (m.sample.index > 0 && m.status == "ok") v.push_back(value_of(m, k));
    KeyStats ks{k, "", statistics(v), failed};
    res.stats.push_back(ks);
  }
  // sensitivity
  for (std::size_t i = 0; i < s.parameters.size(); ++i) {
    for (auto& k : keys) {
      std::vector<double> x, y;
      for (auto& m : res.members)
        if (m.sample.index > 0 && m.status == "ok" && std::isfinite(value_of(m, k))) {
          x.push_back(m.sample.values[i]);
          y.push_back(value_of(m, k));
        }
      double si = first_order_index(x, y);
      res.sobol.push_back({s.parameters[i].path, k, si, std::isfinite(si) && si >= 2.0 / std::sqrt(static_cast<double>(std::max<std::size_t>(1, x.size())))});
      if (s.one_at_a_time) {
        const auto& lo = runs[samples.size() + 2 * i].r, &hi = runs[samples.size() + 2 * i + 1].r;
        res.oat.push_back({s.parameters[i].path, k, value_of(lo, k), value_of(hi, k), value_of(res.members[0], k)});
      }
    }
  }
  // scenarios as factors of the first-order indices (the chosen option)
  for (std::size_t q = 0; q < s.scenarios.size(); ++q)
    for (auto& k : keys) {
      std::vector<double> x, y;
      for (auto& m : res.members)
        if (m.sample.index > 0 && m.status == "ok" && std::isfinite(value_of(m, k))) {
          x.push_back(m.sample.options[q]);
          y.push_back(value_of(m, k));
        }
      double si = first_order_index(x, y);
      res.sobol.push_back({"scenario " + s.scenarios[q].name, k, si, std::isfinite(si) && si >= 2.0 / std::sqrt(static_cast<double>(std::max<std::size_t>(1, x.size())))});
    }
  // probability maps
  std::size_t nok = res.ok();
  for (std::size_t q = 0; q < s.maps.size(); ++q) {
    ProbabilityMap pm;
    pm.spec = s.maps[q];
    const auto& g = map_grids[q];
    pm.x0 = g.x0; pm.y0 = g.y0; pm.dx = g.dx; pm.dy = g.dy; pm.nx = g.nx; pm.ny = g.ny;
    pm.p.assign(g.size(), 0.0);
    for (std::size_t k = 1; k < samples.size(); ++k)
      if (runs[k].r.status == "ok")
        for (std::size_t n = 0; n < g.size(); ++n) pm.p[n] += runs[k].map_hits[q][n];
    for (auto& p : pm.p) p = nok ? p / static_cast<double>(nok) : kNaN;
    res.maps.push_back(std::move(pm));
  }
  if (out_dir.empty()) return res;

  // ---- files ----
  const fs::path e = out_dir / "ensemble";
  fs::create_directories(e / "maps");
  if (bayes) write_posterior_tables(d, post, e);
  const std::string N = std::to_string(s.members), seed = std::to_string(s.seed);
  {
    std::string t = "member,seed,status";
    for (auto& p : s.parameters) t += "," + p.path + " (" + p.unit + ")";
    for (auto& sc : s.scenarios) t += ",scenario " + sc.name;
    for (auto& k : keys) t += "," + csv_text(k);
    t += ",failed_checks,reason,topology\n";
    for (auto& m : res.members) {
      t += (m.sample.index == 0 ? std::string("base") : std::to_string(m.sample.index)) + "," + (m.sample.index == 0 ? std::string("-") : std::to_string(m.sample.seed)) +
           "," + m.status;
      for (double v : m.sample.values) t += "," + fmt_fixed(v, 6);
      for (std::size_t q = 0; q < s.scenarios.size(); ++q) t += "," + s.scenarios[q].options[static_cast<std::size_t>(m.sample.options[q])].name;
      for (auto& k : keys) {
        double v = value_of(m, k);
        t += "," + (std::isfinite(v) ? fmt_shortest(v) : std::string());
      }
      std::string fc;
      for (auto& c : m.failed_checks) fc += (fc.empty() ? "" : "; ") + c;
      t += "," + csv_text(fc) + "," + csv_text(m.reason) + "," + csv_text(m.sample.index == 0 ? std::string("base") : m.topology) + "\n";
    }
    write_file(e / "members.csv", t);
  }
  {
    std::string t = "key,n,failed,mean,sd,p10,p50,p90,min,max,members,seed\n";
    for (auto& k : res.stats)
      t += csv_text(k.key) + "," + std::to_string(k.stats.n) + "," + std::to_string(k.failed) + "," + fmt_fixed(k.stats.mean, 6) + "," + fmt_fixed(k.stats.sd, 6) + "," +
           fmt_fixed(k.stats.p10, 6) + "," + fmt_fixed(k.stats.p50, 6) + "," + fmt_fixed(k.stats.p90, 6) + "," + fmt_fixed(k.stats.min, 6) + "," +
           fmt_fixed(k.stats.max, 6) + "," + N + "," + seed + "\n";
    write_file(e / "statistics.csv", t);
  }
  {
    std::string t = "parameter,key,first_order_index,distinguishable_from_0,members_ok,seed\n";
    for (auto& r : res.sobol)
      t += r.parameter + "," + csv_text(r.key) + "," + fmt_fixed(r.index, 6) + "," + (r.distinguishable ? "yes" : "no") + "," + std::to_string(nok) + "," + seed + "\n";
    write_file(e / "sensitivity_sobol.csv", t);
  }
  if (s.one_at_a_time) {
    std::string t = "parameter,key,value_at_p10,value_at_p90,base,swing\n";
    for (auto& r : res.oat)
      t += r.parameter + "," + csv_text(r.key) + "," + fmt_shortest(r.low) + "," + fmt_shortest(r.high) + "," + fmt_shortest(r.base) + "," +
           fmt_shortest(std::fabs(r.high - r.low)) + "\n";
    write_file(e / "sensitivity_oat.csv", t);
  }
  for (auto& pm : res.maps) {
    Grid2D g(GridSpec{pm.x0, pm.y0, pm.dx, pm.dy, pm.nx, pm.ny});
    g.v = pm.p;
    const std::string nm = pm.spec.kind == "unit" ? "P_" + io::petrel_name(pm.spec.unit) + "_at_" + fmt_fixed(pm.spec.depth, 0) + "m"
                                                  : "P_closure_" + io::petrel_name(pm.spec.horizon);
    write_file(e / "maps" / (nm + ".irap"), io::write_irap_ascii(g, 6));
  }
  std::string readme = "Stratonix ensemble\n\nModel: " + d.model.name + " v" + d.model.version + " (" + d.model.status + ")\n" +
                       "Members: " + N + "; seed: " + seed + "; generator: xoshiro256** seeded per member by SplitMix64(seed, member)\n" +
                       "Members ok: " + std::to_string(nok) + "; failed (refused, build failed or checks failed): " + std::to_string(failed) + " - every one listed in members.csv\n" +
                       "Members whose topology (contact graph) differs from the base case: " + std::to_string(res.topology_differs) + " - the differing contacts in members.csv\n\n" +
                       "members.csv          every member: seed, sampled values, scenario options, status and reason, key numbers\n" +
                       (bayes ? std::string("posterior_*.csv      posterior of the Bayesian conditioning: retained samples, priors, diagnostics (members = these samples; base = least-squares fit)\n") : std::string()) +
                       "statistics.csv       per key number over the ok members: n, mean, sd, P10, P50, P90, min, max (with N and seed)\n" +
                       "sensitivity_sobol.csv first-order variance-based indices from the ensemble (given data, sqrt(N) bins)\n" +
                       (s.one_at_a_time ? std::string("sensitivity_oat.csv  one at a time: each parameter at its P10 and P90, others at base\n") : std::string()) +
                       "maps/*.irap          probability maps (fraction of the ok members)\n" +
                       (o.figures ? std::string("figures/             histograms, first-order indices, tornado charts and probability maps (PNG, SVG, TIFF, PPTX)\n") : std::string()) +
                       (o.base_build ? std::string("../base/             the base case (the description as written): full build, checks, exports and figures\n") : std::string());
  if (o.figures) {
    auto m0 = build_only(dbase, BuildOptions{}).model;
    std::vector<std::pair<std::string, std::string>> kk = {{"Members ok / failed", std::to_string(nok) + " / " + std::to_string(failed)},
                                                           {"Seed", seed}};
    auto c = figure_context(d, *m0, kk);
    const fig::EnsembleLabel lab{nok, failed, s.seed};
    const std::string ev = std::string("Evidence basis: ") + basis_label(m0->info().basis) + ". Ensemble of " + N + " members, seed " + seed + ".";
    auto put = [&](const fig::Figure& f, const std::string& kind, const std::string& name, const std::string& t2, const std::string& cap) {
      fig::write_figure(f, e / "figures", fig::figure_basename(c, kind, name), c.model_name + " - " + t2, cap, ev);
    };
    for (auto& ks : res.stats) {
      std::vector<double> v;
      for (auto& m : res.members)
        if (m.sample.index > 0 && m.status == "ok") v.push_back(value_of(m, ks.key));
      const std::string nm = io::petrel_name(ks.key);
      put(fig::histogram_figure(c, ks.key, v, ks.stats.p10, ks.stats.p50, ks.stats.p90, lab), "histogram", nm, "distribution of " + ks.key,
          "Histogram of " + ks.key + " over " + std::to_string(nok) + " members (seed " + seed + ").");
      std::vector<std::pair<std::string, double>> idx;
      for (auto& r : res.sobol)
        if (r.key == ks.key) idx.push_back({r.parameter, r.index});
      if (!idx.empty())
        put(fig::sobol_figure(c, ks.key, idx, 2.0 / std::sqrt(static_cast<double>(std::max<std::size_t>(1, nok))), lab), "sensitivity", nm,
            "first-order indices of " + ks.key, "First-order variance-based indices of " + ks.key + ".");
      std::vector<fig::TornadoBar> bars;
      double base = kNaN;
      for (auto& r : res.oat)
        if (r.key == ks.key) { bars.push_back({r.parameter, r.low, r.high}); base = r.base; }
      if (!bars.empty())
        put(fig::tornado_figure(c, ks.key, base, bars, lab), "tornado", nm, "sensitivity of " + ks.key, "One-at-a-time sensitivity of " + ks.key + ".");
    }
    for (auto& pm : res.maps) {
      Grid2D g(GridSpec{pm.x0, pm.y0, pm.dx, pm.dy, pm.nx, pm.ny});
      g.v = pm.p;
      const std::string what = pm.spec.kind == "unit" ? pm.spec.unit + " at " + fmt_fixed(pm.spec.depth, 0) + " m depth" : "inside the closure of " + pm.spec.horizon;
      put(fig::probability_map_figure(c, what, g, lab), "probability", io::petrel_name(pm.spec.kind == "unit" ? pm.spec.unit : pm.spec.horizon), "probability map",
          "Probability (fraction of " + std::to_string(nok) + " members) of " + what + ".");
    }
  }
  if (o.base_build) {
    BuildOptions bo;
    bo.figures = o.figures;
    run_build(description, out_dir / "base", bo, "ensemble base case");
  }
  write_file(e / "README.txt", readme);
  Manifest man = Manifest::from_description(d, "ensemble");
  man.extra.push_back({"ensemble.members", N});
  man.extra.push_back({"ensemble.seed", seed});
  man.extra.push_back({"ensemble.generator", "xoshiro256** / SplitMix64"});
  man.add_outputs(out_dir);
  man.write(out_dir / "manifest.json");
  return res;
}

}  // namespace stratonix::unc
