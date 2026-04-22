#pragma once

#include "types.h"
#include <string>
#include <string_view>
#include <ostream>

namespace seq {

class Sequence {
public:
  Sequence() = default;
  explicit Sequence(std::string str);
  Sequence(const char* str);

  std::size_t size() const noexcept;
  bool empty() const noexcept;

  Symbol  operator[](std::size_t i) const;

  Alphabet alphabet() const;

  const std::string& str() const noexcept;

  auto begin() const noexcept { return data_.begin(); }
  auto end()   const noexcept { return data_.end(); }

  bool operator==(const Sequence& other) const = default;

  friend std::ostream& operator<<(std::ostream& os, const Sequence& s);

private:
  std::string data_;
};

} // namespace seq
