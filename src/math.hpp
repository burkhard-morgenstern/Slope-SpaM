#include <cmath>
#include <iterator>
#include <numeric>
#include <optional>

template <class Range, class Proj>
auto mean(Range&& rng, Proj proj = [](auto&& e) { return e; }) {
  using ResultType = std::decay_t<decltype(proj(*std::begin(rng)))>;
  auto sum = ResultType{0};
  for (auto& e : rng) {
    sum += proj(e);
  }
  return ResultType{1} / rng.size() * sum;
}

template <class Range>
auto slope(Range&& rng) {
  using ResultType = std::decay_t<decltype(std::get<1>(*std::begin(rng)))>;
  auto const x_mean = mean(rng, [](auto& e) { return e.first; });
  auto const y_mean = mean(rng, [](auto& e) { return e.second; });
  auto num = ResultType{0};
  auto denom = ResultType{0};
  for (auto&& p : rng) {
    auto const& [x, y] = p;
    num += (x - x_mean) * (y - y_mean);
    denom += (x - x_mean) * (x - x_mean);
  }
  return num / denom;
}
