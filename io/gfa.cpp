#include "gfa.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

namespace seq {

SequenceGraph
read_gfa(std::istream& is)
{
  std::map<VertexId, Symbol> segments;
  std::vector<std::pair<VertexId, VertexId>> links;

  std::string line;
  while (std::getline(is, line)) {
    if (line.empty() || line[0] == '#' || line[0] == 'H')
      continue;

    if (line[0] == 'S') {
      std::istringstream ss(line);
      char type;
      VertexId id;
      std::string label;
      ss >> type >> id >> label;
      if (!label.empty())
        segments[id] = label[0];
    }
    else if (line[0] == 'L') {
      std::istringstream ss(line);
      char type;
      VertexId u, v;
      std::string orient_u, orient_v, overlap;
      // Orientacao (+/-) e overlap (CIGAR) sao lidos e descartados
      // intencionalmente: o SGAP trata o grafo como direcionado simples.
      ss >> type >> u >> orient_u >> v >> orient_v >> overlap;
      links.emplace_back(u, v);
    }
  }

  if (segments.empty())
    return SequenceGraph(0);

  std::size_t vertex_count = segments.rbegin()->first + 1;
  SequenceGraph g(vertex_count);

  for (auto& [id, label] : segments)
    g.set_label(id, label);

  for (auto& [u, v] : links)
    g.add_arc(u, v);

  return g;
}

SequenceGraph
read_gfa(const std::filesystem::path& path)
{
  std::ifstream ifs(path);
  if (!ifs)
    throw std::runtime_error("cannot open file: " + path.string());
  return read_gfa(ifs);
}

} // namespace seq
