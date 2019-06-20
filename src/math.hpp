#include <cmath>
#include <iterator>
#include <numeric>
#include <optional>

#include <range/v3/view.hpp>
#include <range/v3/range/traits.hpp>

template<class Range>
auto mean(Range&& rng)
{
	return ranges::range_value_t<Range>{1} / ranges::size(rng) *
		ranges::accumulate(rng, ranges::range_value_t<Range>{0});
}

template<class Range>
auto slope(Range&& rng)
{
	using ResultType = std::decay_t<decltype(std::get<1>(*ranges::begin(rng)))>;
	namespace rv = ranges::view;
	auto const x_mean = mean(rng | rv::keys);
	auto const y_mean = mean(rng | rv::values);
	auto num = ResultType{0};
	auto denom = ResultType{0};
	for (auto&& p : rng) {
		auto const& [x, y] = p;
		num += (x - x_mean) * (y - y_mean);
		denom += (x - x_mean) * (x - x_mean);
	}
	return static_cast<ResultType>(num / denom);
}
