/// @file linalg.hpp
/// @brief Dense linear solvers for interpolation systems (Eigen is not available in every build;
/// these are small, dependency-free and deterministic).
#pragma once
#include <cstddef>
#include <vector>

namespace stratonix {

/// Solve A x = b for a dense n x n row-major matrix by Gaussian elimination with partial
/// pivoting. Reference implementation: plain loops, the specification for the optimised solver.
/// A and b are overwritten; the solution is returned in b. Throws Error(Build) if A is singular.
void solve_dense_reference(std::vector<double>& A, std::vector<double>& b, std::size_t n);

/// Same result by a blocked right-looking LU factorisation with partial pivoting whose trailing
/// updates run in parallel (fixed chunks, so the result does not depend on the thread count).
/// Pivot choices equal the reference's up to rounding; parity is tested.
void solve_dense_blocked(std::vector<double>& A, std::vector<double>& b, std::size_t n, std::size_t block = 64);

/// Solve with the optimised solver when optimised_kernels() is on and n >= threshold, otherwise
/// with the reference solver.
void solve_dense(std::vector<double>& A, std::vector<double>& b, std::size_t n, std::size_t threshold = 256);

/// Solve A X = B for m right-hand sides (B row-major n x m), by the reference elimination with
/// partial pivoting. A and B are overwritten; the solutions are returned in B. Used where
/// entries of the inverse are needed (leave-one-out errors, stage 3A).
void solve_dense_multi(std::vector<double>& A, std::vector<double>& B, std::size_t n, std::size_t m);

/// Eigen-decomposition of a symmetric 3 x 3 matrix by Jacobi rotations: eigenvalues ascending,
/// eigenvectors as columns (v[3*row + col]).
void symmetric_eigen3(const double a[9], double values[3], double vectors[9]);

}  // namespace stratonix
