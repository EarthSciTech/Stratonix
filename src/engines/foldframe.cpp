// Stage 4A (STX-TN-14 section 2): fold frames. Interface only until Prompt 44.
#include "stratonix/engines/foldframe.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

double FoldFrame::frame_coordinate(const Vec3&) const { throw Error(ErrorCode::NotImplemented, "stage 4A: fold frames (Prompt 44)"); }
double FoldFrame::rotation_angle(double) const { throw Error(ErrorCode::NotImplemented, "stage 4A: fold frames (Prompt 44)"); }
Vec3 FoldFrame::bedding_normal(const Vec3&) const { throw Error(ErrorCode::NotImplemented, "stage 4A: fold frames (Prompt 44)"); }
FoldFrame fit_fold_frame(const std::vector<OrientationDatum>&, const std::vector<OrientationDatum>&, const FoldFrameOptions&) { throw Error(ErrorCode::NotImplemented, "stage 4A: fold frames (Prompt 44)"); }

}  // namespace stratonix
