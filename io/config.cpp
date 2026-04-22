#include "config.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace seq {

static std::string trim(const std::string& s) {
  auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

static std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

static bool parse_bool(const std::string& val) {
  auto v = to_lower(val);
  return v == "true" || v == "1" || v == "yes";
}


ExperimentConfig read_config(const std::string& path) {
  std::ifstream file(path);
  if (!file.is_open())
    throw std::runtime_error("nao foi possivel abrir config: " + path);

  ExperimentConfig cfg;
  std::string line;

  while (std::getline(file, line)) {
    // Remover comentario inline
    auto hash = line.find('#');
    if (hash != std::string::npos)
      line = line.substr(0, hash);

    line = trim(line);
    if (line.empty()) continue;

    auto eq = line.find('=');
    if (eq == std::string::npos) continue;

    auto key = trim(line.substr(0, eq));
    auto val = trim(line.substr(eq + 1));

    if      (key == "algorithm")           cfg.algorithm = to_lower(val);
    else if (key == "jain_method")        cfg.jain_method = to_lower(val);
    else if (key == "greedy_restarts")    cfg.greedy_restarts = std::stoi(val);
    else if (key == "greedy_selection")   cfg.greedy_selection = to_lower(val);
    else if (key == "greedy_seed")        cfg.greedy_seed = static_cast<unsigned>(std::stoul(val));
    else if (key == "beam_width")         cfg.beam_width = std::stoi(val);
    else if (key == "beam_width_mode")   cfg.beam_width_mode = to_lower(val);
    else if (key == "beam_potentials")      cfg.beam_potentials = to_lower(val);
    else if (key == "cost_substitution")  cfg.costs.substitution = std::stoi(val);
    else if (key == "cost_insertion")     cfg.costs.insertion = std::stoi(val);
    else if (key == "cost_deletion")      cfg.costs.deletion = std::stoi(val);
    else if (key == "pli_verbose")        cfg.pli_verbose = parse_bool(val);
    else if (key == "pli_time_limit_s")   cfg.pli_time_limit_s = std::stod(val);
  }

  return cfg;
}

} // namespace seq
