#pragma once

#include <filesystem>
#include <iosfwd>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace spam {

struct unassembled_sequence {
	std::string name;
	std::vector<std::string> reads;
	long double error_rate;
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

	auto error_rate() const
		-> long double;
};

auto background_match_probability(
    spam::sequence const& sequence1,
    spam::sequence const& sequence2)
    -> double;

auto load_directory(std::filesystem::path const& path)
	-> std::optional<std::vector<sequence>>;

auto load_fasta_file(std::filesystem::path const& filename)
	-> std::optional<sequence>;

auto load_multi_fasta_file(std::filesystem::path const& filename)
	-> std::optional<std::vector<assembled_sequence>>;

} // namespace spam
