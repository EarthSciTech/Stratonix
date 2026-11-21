// Stage 4A (STX-TN-14 section 6): contact graph and watertightness. Interface only until Prompt 46.
#include "stratonix/repr/topology.hpp"

#include "stratonix/core/error.hpp"

namespace stratonix {

std::string ContactGraph::json() const { throw Error(ErrorCode::NotImplemented, "stage 4A: topology (Prompt 46)"); }
ContactGraph contact_graph(const Model&, double) { throw Error(ErrorCode::NotImplemented, "stage 4A: topology (Prompt 46)"); }
bool same_topology(const ContactGraph&, const ContactGraph&, std::vector<std::string>*) { throw Error(ErrorCode::NotImplemented, "stage 4A: topology (Prompt 46)"); }
WatertightReport watertightness(const Model&, double, double) { throw Error(ErrorCode::NotImplemented, "stage 4A: topology (Prompt 46)"); }

}  // namespace stratonix
