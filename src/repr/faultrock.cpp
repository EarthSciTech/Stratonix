// Stage 2F (STX-TN-10). Implemented in Prompts 32 and 33.
#include "stratonix/repr/faultrock.hpp"

namespace stratonix::repr {

FaultRockPoint fault_rock_at(const Model&, int, const Vec3&, const FaultRockInput&, const StressInput&) {
  FaultRockPoint p;
  p.missing = "fault-rock attributes come with stage 2F (Prompt 32)";
  return p;
}

CurvatureGrids curvature(const Model&, int, const GridSpec& g, double) {
  CurvatureGrids c{Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g), Grid2D(g)};
  return c;
}

}  // namespace stratonix::repr
