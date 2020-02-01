#pragma once

#include <string>
#include <vector>

namespace spam {

class pattern {
  std::string bits;
  std::vector<size_t> indices;

 public:
  explicit pattern(std::string bits);

  auto reduce(size_t k) const -> pattern;

  auto begin() const { return indices.begin(); }

  auto end() const { return indices.end(); }

  auto size() const -> size_t { return bits.size(); }

  auto weight() const -> size_t { return indices.size(); }

  auto operator[](size_t n) const -> size_t { return indices[n]; }

  auto operator==(pattern const& rhs) const -> bool {
    return indices == rhs.indices;
  }
};

}  // namespace spam
