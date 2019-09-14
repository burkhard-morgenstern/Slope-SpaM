#include "sequence.hpp"

#include <fstream>

#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

namespace fs = std::filesystem;
namespace rv = ranges::view;

namespace spam {

auto operator>>(std::istream& is, assembled_sequence& seq)
    -> std::istream&
{
    std::getline(is, seq.name);
    std::string nucleotides;
    std::getline(is, nucleotides, '>');
    nucleotides.erase(std::remove_if(nucleotides.begin(), nucleotides.end(),
        [](const char c) { return isspace(c); }), nucleotides.end());
    seq.nucleotides = std::move(nucleotides);
    return is;
}

auto sequence::name() const
    -> std::string
{
    if (std::holds_alternative<assembled_sequence>(*this)) {
        return std::get<assembled_sequence>(*this).name;
    } else {
        return std::get<unassembled_sequence>(*this).name;
    }
}

auto sequence::size() const
    -> size_t
{
    if (std::holds_alternative<assembled_sequence>(*this)) {
        return std::get<assembled_sequence>(*this).size();
    } else {
        return ranges::accumulate(
            std::get<unassembled_sequence>(*this).reads
                | rv::transform(
                    [](auto&& read) {
                        return read.size();
                    }),
            size_t{0});
    }
}

auto sequence::adjusted_size(size_t wordlength) const
    -> size_t
{
    if (std::holds_alternative<assembled_sequence>(*this)) {
        return std::get<assembled_sequence>(*this).size() - wordlength;
    } else {
        auto& seq = std::get<unassembled_sequence>(*this);
        return ranges::accumulate(
            seq.reads
                | rv::transform(
                    [](auto&& read) {
                        return read.size();
                    }),
            size_t{0}) -
            seq.reads.size() * wordlength;
    }
}

auto load_directory(fs::path const& path)
    -> std::optional<std::vector<sequence>>
{
    if (!fs::is_directory(path)) {
        return {};
    }

	auto result = std::vector<sequence>{};
    for (auto& entry : fs::directory_iterator(path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".fasta") {
            result.push_back(*load_fasta_file(entry));
        }
    }

    return {result};
}

auto load_fasta_file(fs::path const& filename)
    -> std::optional<sequence>
{
    auto file = std::ifstream(filename);
    if (!file.is_open()) {
        return {};
    }
	auto reads = std::vector<assembled_sequence>{};
	auto dummy = std::string{};
    std::getline(file, dummy, '>');
    while (!file.eof()) {
		auto seq = assembled_sequence{};
		file >> seq;
		reads.push_back(std::move(seq));
	}

    if (reads.size() == 1) {
        return {sequence{reads[0]}};
    } else {
        auto result = unassembled_sequence{};
        result.name = filename.stem();
        result.reads.resize(reads.size());
        for (auto i = size_t{0}; i < reads.size(); ++i) {
            result.reads[i] = std::move(reads[i].nucleotides);
        }
        return {sequence{result}};
    }
}

auto load_multi_fasta_file(fs::path const& filename)
    -> std::optional<std::vector<assembled_sequence>>
{
    auto file = std::ifstream(filename);
    if (!file.is_open()) {
        return {};
    }
	auto result = std::vector<assembled_sequence>{};
	auto dummy = std::string{};
    std::getline(file, dummy, '>');
    while (!file.eof()) {
		auto seq = assembled_sequence{};
		file >> seq;
		result.push_back(std::move(seq));
	}
	return {result};
}

} // namespace spam
