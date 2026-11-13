/// @file roots.hpp
/// @brief Bracketed scalar root finding (Illinois variant of regula falsi).
#pragma once
#include <cmath>

namespace stratonix {

/// Root of f in [a, b] given f(a) = fa and f(b) = fb of opposite signs. Converges superlinearly
/// and never leaves the bracket. Stops when |f| < ftol or the bracket is shorter than xtol.
template <class F>
double illinois(F&& f, double a, double b, double fa, double fb, double xtol = 1e-9, double ftol = 0,
                int maxit = 200) {
  if (fa == 0) return a;
  if (fb == 0) return b;
  int side = 0;
  double c = a;
  for (int it = 0; it < maxit; ++it) {
    c = (a * fb - b * fa) / (fb - fa);
    if (!(c > std::fmin(a, b) && c < std::fmax(a, b))) c = 0.5 * (a + b);
    double fc = f(c);
    if (fc == 0 || std::fabs(fc) <= ftol || std::fabs(b - a) < xtol) return c;
    if ((fc > 0) == (fa > 0)) {
      a = c; fa = fc;
      if (side == -1) fb *= 0.5;
      side = -1;
    } else {
      b = c; fb = fc;
      if (side == 1) fa *= 0.5;
      side = 1;
    }
  }
  return c;
}

/// Expand a bracket around x0 (step, doubling) until f changes sign, then solve. Returns NaN if
/// no sign change is found within `max_span` of x0.
template <class F>
double expand_and_solve(F&& f, double x0, double step, double max_span, double xtol = 1e-9) {
  double f0 = f(x0);
  if (f0 == 0) return x0;
  double a = x0, fa = f0;
  for (double s = step; s <= max_span; s *= 2) {
    double b = x0 + s, fb = f(b);
    if ((fb > 0) != (f0 > 0)) return illinois(f, a, b, fa, fb, xtol);
    double c = x0 - s, fc = f(c);
    if ((fc > 0) != (f0 > 0)) return illinois(f, c, a, fc, fa, xtol);
  }
  return std::nan("");
}

}  // namespace stratonix
