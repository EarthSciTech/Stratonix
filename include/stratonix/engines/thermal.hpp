// Stage 4B (STX-TN-15 section 6, decision 20): temperature and overpressure along 1D columns.
#pragma once
#include <string>
#include <vector>

#include "stratonix/core/geometry.hpp"

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

// ---- columns through the burial history (models) -------------------------------------------------
class Structure;
/// Per model unit (index): grain conductivity, heat production, and the pressure laws; units
/// without values take the deepest given unit's.
struct ColumnProps {
  std::vector<double> conductivity, heat_production;      ///< W/m/K, W/m3
  std::vector<int> perm_law;                              ///< 0 none, 1 exponential, 2 Kozeny-Carman
  std::vector<double> k0, b, phi0;                        ///< m2, -, -
};
struct ThermalColumnSpec {
  double t_surface = 10;            ///< degC
  double q_base = 0.06;             ///< W/m2 (constant, or the history below)
  std::vector<std::pair<double, double>> q_history;  ///< (age Ma, W/m2), interpolated in time
  double k_water = 0.6;             ///< W/m/K
  double column_depth = 0;          ///< m below the surface (0: twice the deepest horizon, at least 5 km)
  bool vitrinite = false;
  ColumnProps props;
};
/// Temperature and EASY%Ro of each horizon through the steps of the burial history at (x, y).
struct ThermalHistory {
  std::vector<std::string> labels;
  std::vector<double> ages;                     ///< Ma (NaN without ages)
  std::vector<std::vector<double>> temperature; ///< [step][horizon] degC (NaN before deposition)
  std::vector<std::vector<double>> ro;          ///< [step][horizon] %Ro (empty without vitrinite)
  std::vector<double> present_depth, present_temperature;  ///< present profile (m, degC)
  double surface_flow_present = 0;              ///< W/m2, from the present gradient at the top
  double balance_error = 0;                     ///< steady present column: |q_s - (q_b + sum A h)| / q_s
};
ThermalHistory thermal_history(const Structure& st, const ThermalColumnSpec& s, double x, double y);

struct PressureColumnSpec {
  double viscosity = 1e-3, rho_w = 1030, rho_g = 2650, g = 9.81;
  ColumnProps props;
};
/// Excess pore pressure of each horizon through the steps (disequilibrium compaction, Gibson 1958;
/// Bethke 1985), and the present profile with its hydrostatic and lithostatic bounds.
struct PressureHistory {
  std::vector<std::string> labels;
  std::vector<double> ages;
  std::vector<std::vector<double>> overpressure;  ///< [step][horizon] Pa
  std::vector<double> present_depth, present_pore, present_hydrostatic, present_lithostatic;  ///< m, Pa
};
PressureHistory pressure_history(const Structure& st, const PressureColumnSpec& s, double x, double y);

/// EASY%Ro (Sweeney and Burnham 1990) of a time-temperature path (Ma, degC), oldest first.
double easy_ro(const std::vector<std::pair<double, double>>& path);

}  // namespace stratonix
