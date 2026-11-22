/// @file units.hpp
/// @brief Units accepted in model descriptions and their conversion to SI.
#pragma once
#include <optional>
#include <string>
#include <string_view>

namespace stratonix {

/// Physical dimension of a parameter.
enum class Dimension { Length, Angle, Dimensionless, Ordinal, Flag, Area, InverseLength, Time, Pressure, PressureGradient, Permeability,
                       Density, Acceleration, MagneticField,  ///< stage 4A: kg/m3, mGal, nT
                       Temperature, ThermalConductivity, HeatFlow, HeatProduction, Viscosity,  ///< stage 4B: degC, W/m/K, mW/m2, uW/m3, mPa.s
                       Any };

struct UnitDef {
  std::string_view name;
  Dimension dim;
  double to_si;  ///< multiply by this to get SI (m, rad, 1)
};

/// Look up a unit by name ("m", "km", "ft", "deg", "rad", "fraction", "ratio", "percent",
/// "ordinal", "flag"). Returns nullopt if unknown.
std::optional<UnitDef> find_unit(std::string_view name);
const char* dimension_name(Dimension d);

}  // namespace stratonix
