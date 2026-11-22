// Stage 4B (STX-TN-15 section 6, decision 20): temperature and overpressure along 1D columns.
#pragma once
#include <vector>

namespace stratonix {

struct ThermalLayer { double thickness = 0, conductivity = 2, heat_production = 0; };  ///< m, W/m/K, W/m^3
/// Steady temperatures at the top of the column and at each layer base, with surface temperature
/// and surface heat flow (W/m^2, positive upward).
std::vector<double> steady_temperature(const std::vector<ThermalLayer>& layers, double t_surface, double q_surface);

struct TransientSpec {
  double length = 1000, diffusivity = 1e-6;   ///< m, m^2/s
  double t_top = 0, t_base = 100, t_initial = 0;
  double time = 0;                            ///< s after the base step
  int nodes = 201;
};
/// Temperatures at equally spaced nodes (top first) after `time`.
std::vector<double> transient_temperature(const TransientSpec& s);

struct GibsonSpec {
  double rate = 1e-11;          ///< sedimentation rate (m/s)
  double duration = 0;          ///< s
  double permeability = 1e-18;  ///< m^2
  double viscosity = 1e-3;      ///< Pa s
  double compressibility = 1e-8;  ///< 1/Pa (mv)
  double rho_s = 2600, rho_w = 1000, g = 9.81;
  int nodes = 200;
};
/// Overpressure (Pa) at the base of a column deposited at a constant rate (Gibson 1958).
double gibson_overpressure_base(const GibsonSpec& s);

}  // namespace stratonix
