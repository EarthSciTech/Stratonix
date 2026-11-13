// Dense solvers and a 3x3 symmetric eigen-solver.
#include "stratonix/core/linalg.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "stratonix/core/error.hpp"
#include "stratonix/core/util.hpp"

namespace stratonix {

namespace {
void singular(std::size_t k) {
  throw Error(ErrorCode::Build, "interpolation system is singular at column " + std::to_string(k) +
                                    " (duplicate or collinear data, or too few data for the drift terms)");
}
}  // namespace

void solve_dense_reference(std::vector<double>& A, std::vector<double>& b, std::size_t n) {
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t piv = k;
    double best = std::fabs(A[k * n + k]);
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::fabs(A[i * n + k]) > best) { best = std::fabs(A[i * n + k]); piv = i; }
    if (!(best > 0)) singular(k);
    if (piv != k) {
      for (std::size_t j = 0; j < n; ++j) std::swap(A[k * n + j], A[piv * n + j]);
      std::swap(b[k], b[piv]);
    }
    for (std::size_t i = k + 1; i < n; ++i) {
      double f = A[i * n + k] / A[k * n + k];
      A[i * n + k] = f;
      for (std::size_t j = k + 1; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
      b[i] -= f * b[k];
    }
  }
  for (std::size_t k = n; k-- > 0;) {
    double s = b[k];
    for (std::size_t j = k + 1; j < n; ++j) s -= A[k * n + j] * b[j];
    b[k] = s / A[k * n + k];
  }
}

void solve_dense_blocked(std::vector<double>& A, std::vector<double>& b, std::size_t n, std::size_t nb) {
  std::vector<std::size_t> perm(n);
  for (std::size_t i = 0; i < n; ++i) perm[i] = i;
  for (std::size_t k0 = 0; k0 < n; k0 += nb) {
    std::size_t k1 = std::min(n, k0 + nb);
    // Panel factorisation (columns k0..k1-1), rows k0..n-1, with row swaps across the full width.
    for (std::size_t k = k0; k < k1; ++k) {
      std::size_t piv = k;
      double best = std::fabs(A[k * n + k]);
      for (std::size_t i = k + 1; i < n; ++i)
        if (std::fabs(A[i * n + k]) > best) { best = std::fabs(A[i * n + k]); piv = i; }
      if (!(best > 0)) singular(k);
      if (piv != k) {
        for (std::size_t j = 0; j < n; ++j) std::swap(A[k * n + j], A[piv * n + j]);
        std::swap(b[k], b[piv]);
        std::swap(perm[k], perm[piv]);
      }
      double inv = 1.0 / A[k * n + k];
      for (std::size_t i = k + 1; i < n; ++i) {
        double f = A[i * n + k] * inv;
        A[i * n + k] = f;
        for (std::size_t j = k + 1; j < k1; ++j) A[i * n + j] -= f * A[k * n + j];
      }
    }
    if (k1 >= n) break;
    // U12 = L11^-1 A12 (rows k0..k1-1, columns k1..n-1).
    for (std::size_t k = k0; k < k1; ++k)
      for (std::size_t i = k + 1; i < k1; ++i) {
        double f = A[i * n + k];
        if (f != 0)
          for (std::size_t j = k1; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
      }
    // Trailing update A22 -= L21 U12, rows in parallel (each row by one call).
    parallel_for(n - k1, 16, [&](std::size_t r0, std::size_t r1) {
      for (std::size_t r = r0; r < r1; ++r) {
        std::size_t i = k1 + r;
        double* ai = &A[i * n];
        for (std::size_t k = k0; k < k1; ++k) {
          double f = ai[k];
          if (f == 0) continue;
          const double* ak = &A[k * n];
          for (std::size_t j = k1; j < n; ++j) ai[j] -= f * ak[j];
        }
      }
    });
  }
  // Forward substitution (unit L) and back substitution (U).
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= A[i * n + k] * b[k];
    b[i] = s;
  }
  for (std::size_t k = n; k-- > 0;) {
    double s = b[k];
    for (std::size_t j = k + 1; j < n; ++j) s -= A[k * n + j] * b[j];
    b[k] = s / A[k * n + k];
  }
}

void solve_dense(std::vector<double>& A, std::vector<double>& b, std::size_t n, std::size_t threshold) {
  if (optimised_kernels() && n >= threshold) solve_dense_blocked(A, b, n);
  else solve_dense_reference(A, b, n);
}

void symmetric_eigen3(const double a_in[9], double values[3], double vectors[9]) {
  double a[9];
  std::copy(a_in, a_in + 9, a);
  double v[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = std::fabs(a[1]) + std::fabs(a[2]) + std::fabs(a[5]);
    if (off < 1e-300) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        double apq = a[3 * p + q];
        if (std::fabs(apq) < 1e-300) continue;
        double theta = 0.5 * (a[3 * q + q] - a[3 * p + p]) / apq;
        double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 3; ++k) {  // A = J^T A J
          double akp = a[3 * k + p], akq = a[3 * k + q];
          a[3 * k + p] = c * akp - s * akq;
          a[3 * k + q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          double apk = a[3 * p + k], aqk = a[3 * q + k];
          a[3 * p + k] = c * apk - s * aqk;
          a[3 * q + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          double vkp = v[3 * k + p], vkq = v[3 * k + q];
          v[3 * k + p] = c * vkp - s * vkq;
          v[3 * k + q] = s * vkp + c * vkq;
        }
      }
  }
  int idx[3] = {0, 1, 2};
  std::sort(idx, idx + 3, [&](int x, int y) { return a[4 * x] < a[4 * y]; });
  for (int c = 0; c < 3; ++c) {
    values[c] = a[4 * idx[c]];
    for (int r = 0; r < 3; ++r) vectors[3 * r + c] = v[3 * r + idx[c]];
  }
}

}  // namespace stratonix
