#include <algorithm>
#include "string.hpp"

auto split(std::string const& s, char seperator) -> std::vector<std::string> {
  if (s == "") {
    return {};
  }

  auto output = std::vector<std::string>{};
  auto prev_pos = size_t{0};
  auto pos = size_t{0};
  while ((pos = s.find(seperator, pos)) != std::string::npos) {
    auto substring = s.substr(prev_pos, pos - prev_pos);
    output.push_back(substring);
    prev_pos = ++pos;
  }
  output.push_back(s.substr(prev_pos, pos - prev_pos));
  output.erase(std::remove_if(output.begin(), output.end(),
                              [](std::string const& str) {
                                return str == "";
                              }),
               output.end());
  return output;
}
