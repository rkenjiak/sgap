#pragma once

#include "../core/sequence.h"
#include <filesystem>
#include <iostream>

namespace seq {

Sequence read_fastq(const std::filesystem::path& path);
Sequence read_fastq(std::istream& is);

} // namespace seq
