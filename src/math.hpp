#include <cmath>
#include <iterator>
#include <numeric>
#include <optional>

template<class InputIt, class ResultType = typename InputIt::value_type>
auto mean(InputIt first, InputIt last)
	-> ResultType
{
	return ResultType{1} / std::distance(first, last) *
		std::accumulate(first, last, ResultType{0});
}

template<class ForwardIt1, class ForwardIt2,
	class ResultType = std::common_type_t<
		typename ForwardIt1::value_type,
		typename ForwardIt2::value_type>>
auto slope(ForwardIt1 x_first, ForwardIt1 x_last,
	ForwardIt2 y_first, ForwardIt2 y_last)
	-> std::optional<ResultType>
{
	if (std::distance(x_first, x_last) != std::distance(y_first, y_last)) {
		return {};
	}
	auto const x_mean = mean(x_first, x_last);
	auto const y_mean = mean(y_first, y_last);
	auto num = ResultType{0};
	auto denom = ResultType{0};
	auto x_it = x_first;
	auto y_it = y_first;
	for (; x_it != x_last; ++x_it, ++y_it){
		num += (*x_it - x_mean) * (*y_it - y_mean);
		denom += (*x_it - x_mean) * (*x_it - x_mean);
	}
	return {num / denom};
}
