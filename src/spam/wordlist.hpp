#pragma once

#include "pattern.hpp"
#include "sequence.hpp"

#include <range/v3/view/transform.hpp>

namespace spam {

using word_t = __uint128_t;

class wordlist {
	std::vector<word_t> words;
	spam::pattern pattern;

public:
	static auto max_wordsize()
		-> size_t;

	static auto kmin(sequence const&)
		-> size_t;
	static auto kmin(std::vector<sequence> const&)
		-> size_t;

	static auto kmax(sequence const&)
		-> size_t;
	static auto kmax(std::vector<sequence> const&)
		-> size_t;

public:
	wordlist(
		spam::sequence const& sequence,
		spam::pattern pattern);

	auto size() const
	{
		return words.size();
	}

	auto begin() const
	{
		return words.begin();
	}

	auto end() const
	{
		return words.end();
	}

	auto reduce(size_t k) const
	{
		auto reduction_pattern = std::numeric_limits<word_t>::max();
		reduction_pattern <<= 8 * sizeof(word_t) - 2 * k;
		return *this
			| ranges::view::transform(
				[reduction_pattern = reduction_pattern](auto word) {
					return word & reduction_pattern;
				});
	}
};

} // namespace spam
