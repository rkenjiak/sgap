#pragma once

#include "../core/sequence_graph.h"
#include <filesystem>
#include <iostream>
#include <string>

namespace seq {

SequenceGraph read_gfa(const std::filesystem::path& path);
SequenceGraph read_gfa(std::istream& is);

} // namespace seq
