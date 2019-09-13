#pragma once

#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <range/v3/view/transform.hpp>

#include <threadpool/ThreadPool.h>

namespace spam {

class pattern {
	std::string bits;
	std::vector<size_t> indices;

public:
	pattern(std::string bits);

	auto reduce(size_t k) const
		-> pattern;

	auto begin() const
	{
		return indices.begin();
	}

	auto end() const
	{
		return indices.end();
	}

	auto size() const
		-> size_t
	{
		return bits.size();
	}

	auto weight() const
		-> size_t
	{
		return indices.size();
	}

	auto operator[](size_t n) const
		-> size_t
	{
		return indices[n];
	}

	auto operator==(pattern const& rhs) const
		-> bool
	{
		return indices == rhs.indices;
	}
};

struct unassembled_sequence {
	std::string name;
	std::vector<std::string> reads;
};

struct assembled_sequence {
	std::string name;
	std::string nucleotides;

	auto size() const
		-> size_t
	{
		return nucleotides.size();
	}

	auto begin() const
	{
		return nucleotides.begin();
	}

	auto end() const
	{
		return nucleotides.end();
	}
};

auto operator>>(std::istream& is, assembled_sequence& seq)
	-> std::istream&;

struct sequence
	: public std::variant<unassembled_sequence, assembled_sequence>
{
	sequence(unassembled_sequence seq)
		: std::variant<unassembled_sequence, assembled_sequence>{seq}
	{};
	sequence(assembled_sequence seq)
		: std::variant<unassembled_sequence, assembled_sequence>{seq}
	{};

	auto name() const
		-> std::string;

	auto size() const
		-> size_t;

	auto adjusted_size(size_t wordlength) const
		-> size_t;
};

auto load_directory(std::filesystem::path const& path)
	-> std::optional<std::vector<sequence>>;

auto load_fasta_file(std::filesystem::path const& filename)
	-> std::optional<sequence>;

auto load_multi_fasta_file(std::filesystem::path const& filename)
	-> std::optional<std::vector<assembled_sequence>>;

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

class distance_matrix {
	std::vector<spam::sequence> sequences;
	spam::pattern pattern;
	std::vector<size_t> wordlengths;

	std::vector<wordlist> wordlists;
	std::vector<std::vector<double>> matrix;

	std::optional<std::shared_ptr<ThreadPool>> threadpool;

public:
	distance_matrix(
		std::vector<spam::sequence> sequences,
		spam::pattern pattern,
		std::vector<size_t> wordlengths,
		std::optional<std::shared_ptr<ThreadPool>> threadpool = {});

	auto size() const
		-> size_t;

	auto column(size_t i) const
		-> std::pair<spam::sequence const&, std::vector<double> const&>;

private:
	void calculate();

	void initialize_matrix();
	void create_wordlists();
	void create_wordlists_seq();
	void create_wordlists_par();
	void calculate_matrix();
	void calculate_matrix_seq();
	void calculate_matrix_par();

	auto calculate_element(size_t i, size_t j) const
		-> std::pair<double, double>;
};

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix);

} // namespace spam
