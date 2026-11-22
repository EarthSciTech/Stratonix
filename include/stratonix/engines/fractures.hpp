// Stage 4B (STX-TN-15 section 5, decision 19): fracture intensity, discrete fracture networks,
// mechanical stratigraphy and laterally varying stress.
#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "stratonix/core/geometry.hpp"

namespace stratonix {

struct FractureLaw { double p32_0 = 0, a_k = 0, a_e = 0, p32_max = 10; };
/// P32 (1/m) from the maximum absolute principal curvature (1/m) and maximum strain.
double fracture_intensity(double k_max, double e_max, const FractureLaw& law);

struct FractureSet { Vec3 mean_pole{0, 0, 1}; double kappa = 20, exponent = 2.5, rmin = 1, rmax = 50, weight = 1; };
struct Disk { Vec3 c, n; double r = 0; int set = 0; };
/// Baecher disks in a box with expected P32; refused (empty, reason) above `limit` disks.
struct DfnResult { bool ok = false; std::string reason; std::vector<Disk> disks; };
DfnResult sample_dfn(const BBox& box, double p32, const std::vector<FractureSet>& sets, std::uint64_t seed, std::size_t limit = 1000000);
/// Realised P32 of disks (area of disk parts inside the box / box volume).
double realised_p32(const std::vector<Disk>& d, const BBox& box);

/// Horizontal stresses of a layer (Thiercelin and Plumb 1994): returns (sigma_h, sigma_H) in Pa.
std::pair<double, double> horizontal_stress(double sigma_v, double pore, double alpha, double young, double poisson, double e_h, double e_H);

}  // namespace stratonix
