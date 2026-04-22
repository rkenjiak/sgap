#pragma once

#include "../core/sequence.h"
#include <vector>
#include <numeric>
#include <algorithm>

namespace seq {

// Distancia de edicao (Levenshtein) com custos unitarios.
// Complexidade: O(n*m) tempo, O(min(n,m)) espaco.
inline int
edit_distance(const Sequence& a, const Sequence& b)
{
  const std::size_t n = a.size();
  const std::size_t m = b.size();

  // Usa a sequencia menor como coluna para minimizar memoria
  if (n < m) return edit_distance(b, a);

  std::vector<int> dp(m + 1);
  std::iota(dp.begin(), dp.end(), 0); // dp[j] = j

  for (std::size_t i = 1; i <= n; ++i) {
    int prev = dp[0]; // dp[i-1][j-1]
    dp[0] = static_cast<int>(i);

    for (std::size_t j = 1; j <= m; ++j) {
      int temp = dp[j];
      if (a[i - 1] == b[j - 1]) {
        dp[j] = prev;
      } else {
        dp[j] = 1 + std::min({prev, dp[j], dp[j - 1]});
      }
      prev = temp;
    }
  }

  return dp[m];
}

} // namespace seq
