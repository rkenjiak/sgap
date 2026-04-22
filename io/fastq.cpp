#include "fastq.h"

#include <fstream>
#include <stdexcept>

namespace seq {

// Premissa: 1 read por arquivo. O parser le o primeiro bloco de sequencia
// apos o header '@' e para no separador '+'; reads subsequentes sao ignorados.
Sequence
read_fastq(std::istream& is)
{
  std::string line;
  std::string data;

  // Find header line (starts with '@')
  while (std::getline(is, line)) {
    if (!line.empty() && line[0] == '@')
      break;
  }

  // Read sequence lines until '+' separator or next header
  while (std::getline(is, line)) {
    if (line.empty())
      continue;
    if (line[0] == '+' || line[0] == '@')
      break;
    data += line;
  }

  return Sequence(std::move(data));
}

Sequence
read_fastq(const std::filesystem::path& path)
{
  std::ifstream ifs(path);
  if (!ifs)
    throw std::runtime_error("cannot open file: " + path.string());
  return read_fastq(ifs);
}

} // namespace seq
