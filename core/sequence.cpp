#include "sequence.h"

namespace seq {

Sequence::Sequence(std::string str)
  : data_(std::move(str)) {}

Sequence::Sequence(const char* str)
  : data_(str) {}

std::size_t
Sequence::size() const noexcept
{
  return data_.size();
}

bool
Sequence::empty() const noexcept
{
  return data_.empty();
}

Symbol
Sequence::operator[](std::size_t i) const
{
  return data_[i];
}

Alphabet
Sequence::alphabet() const
{
  return Alphabet(data_.begin(), data_.end());
}

const std::string&
Sequence::str() const noexcept
{
  return data_;
}

std::ostream&
operator<<(std::ostream& os, const Sequence& s)
{
  return os << s.data_;
}

} // namespace seq
