#include "gurobi_sol.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace seq {

std::unordered_map<std::string, double>
read_gurobi_sol(const std::string& path)
{
  std::ifstream in(path);
  if (!in.is_open())
    throw std::runtime_error("nao foi possivel abrir .sol: " + path);

  std::unordered_map<std::string, double> values;
  std::string line;

  while (std::getline(in, line)) {
    auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) continue;
    if (line[first] == '#') continue;

    std::istringstream iss(line);
    std::string name;
    double      val;
    if (iss >> name >> val)
      values.emplace(std::move(name), val);
  }

  return values;
}

void
write_gurobi_sol(
    const std::string& path,
    const std::vector<std::pair<std::string, double>>& assignments,
    const std::vector<std::string>& header_comments)
{
  auto parent = std::filesystem::path(path).parent_path();
  if (!parent.empty())
    std::filesystem::create_directories(parent);

  std::ofstream out(path);
  if (!out.is_open())
    throw std::runtime_error("nao foi possivel escrever .sol: " + path);

  for (const auto& c : header_comments)
    out << "# " << c << "\n";

  for (const auto& [name, value] : assignments)
    out << name << " " << value << "\n";
}

} // namespace seq
