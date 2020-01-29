#include "pattern.hpp"

#include <stdexcept>

namespace spam {

pattern::pattern(std::string _bits) : bits(std::move(_bits)) {
  auto end = static_cast<ptrdiff_t>(bits.size());
  for (; end >= 0 && bits[end] != '1'; --end)
    ;
  bits.resize(end + 1);
  for (size_t i = 0; i < bits.size(); ++i) {
    if (bits[i] == '1') {
      indices.push_back(i);
    }
  }
}

auto pattern::reduce(size_t k) const -> pattern {
  if (weight() < k) {
    throw std::invalid_argument{"weight < k"};
  }

  if (weight() == k) {
    return *this;
  }

  return pattern{{bits.begin(), bits.begin() + indices[k]}};
}

}  // namespace spam
