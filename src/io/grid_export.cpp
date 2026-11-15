// Grid exporters (stage 3B, STX-TN-12 section 6): GRDECL here; RESQML through FESAPI where it is
// available. Opened by Prompt 37, implemented by Prompt 38.
#include <stdexcept>

#include "stratonix/io/export.hpp"

namespace stratonix::io {

std::string write_grdecl(const repr::CornerPointGrid&, const repr::PropertySet*, const repr::TransmissibilityResult*,
                         const GrdeclInfo&, int) {
  throw std::runtime_error("not implemented yet: GRDECL exporter (STX-TN-12, Prompt 38)");
}

}  // namespace stratonix::io
