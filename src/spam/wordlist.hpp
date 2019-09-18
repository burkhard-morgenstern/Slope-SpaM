#pragma once

#include "pattern.hpp"
#include "sequence.hpp"

namespace spam {

using word_t = uint64_t;

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
};

auto calculate_matches(
    spam::wordlist const& wordlist1,
    spam::wordlist const& wordlist2,
    size_t k)
    -> size_t;

} // namespace spam
