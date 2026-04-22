#pragma once

#include "solver_interface.h"
#include "../io/config.h"
#include <optional>

namespace seq {

std::optional<SolverResult>
solve_pli(const AlignmentGraph& ag, const ExperimentConfig& cfg);

} // namespace seq
