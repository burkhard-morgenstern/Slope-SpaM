#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include <fmt/format.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/iterator.hpp>
#include <range/v3/view.hpp>

#include "config.hpp"
#include "spam/distance_matrix.hpp"

namespace fs = std::filesystem;
namespace rv = ranges::view;

class application {
	spam::config config;
	std::shared_ptr<ThreadPool> threadpool;

public:
	application(spam::config config)
		: config(std::move(config)),
		threadpool(std::make_shared<ThreadPool>(
			std::thread::hardware_concurrency()))
	{}

private:
	auto output_path(fs::path const& input_path)
		-> fs::path
	{
		return !fs::is_directory(input_path)
			? fs::path{input_path}.replace_extension(".dmat")
			: fs::path{input_path.string() + ".dir.dmat"};
	}

	auto input_output_mapping()
	{
		return config.in
			| rv::transform(
				[&](auto&& input_path) {
					return std::make_pair(
						input_path,
						(config.out == "" || config.in.size() != 1)
							? output_path(input_path)
							: config.out);
				});
	}

public:
	auto exec()
		-> int
	{
		if (!config.multi_fasta_as_reads) {
			for (auto const& [input_path, output_path] : input_output_mapping()) {
				fmt::print("Processing {} -> {}{}",
					input_path.string(), output_path.string(), config.in.size() < 100 ? "\n" : "\r");
				auto sequences = load_sequences(input_path);
				if (sequences.size() == 0) {
					fmt::print(stderr, "Empty input path \"{}\".\n", input_path.string());
				} else {
					auto file = std::ofstream{output_path};
					file << process(sequences);
					file.close();
				}
			}
		} else {
			fmt::print("Processing ");
			ranges::copy(config.in, ranges::ostream_joiner(std::cout, ' '));
			fmt::print(" -> {}\n", config.out.string());
			auto sequences = config.in
				| rv::transform(
					[&](auto&& input_path) {
						return *spam::load_fasta_file(input_path);
					})
				| ranges::to<std::vector<spam::sequence>>();
			auto file = std::ofstream{config.out};
			file << process(sequences);
			file.close();
		}

		return 0;
	}

private:
	auto load_sequences(fs::path const& path)
		-> std::vector<spam::sequence>
	{
		if (fs::is_directory(path)) {
			return *spam::load_directory(path);
		} else {
			if (config.multi_fasta_as_reads) {
				return std::vector<spam::sequence>{*spam::load_fasta_file(path)};
			} else {
				auto assembled_sequences = *spam::load_multi_fasta_file(path);
				return assembled_sequences
					| rv::transform(
						[](auto&& seq) {
							return spam::sequence{std::move(seq)};
						})
					| ranges::to<std::vector<spam::sequence>>();
			}
		}
	}

	auto process(std::vector<spam::sequence> const& sequences)
		-> spam::distance_matrix
	{
		auto wordlengths = config.wordlengths;
		if (wordlengths.empty()) {
			wordlengths.push_back(spam::wordlist::kmin(sequences));
			wordlengths.push_back(spam::wordlist::kmax(sequences));
		}
		size_t kmax = *ranges::max_element(wordlengths);
		auto pattern = config.pattern;
		return spam::distance_matrix(
			std::move(sequences),
			pattern,
			wordlengths,
			threadpool);
	}
};

int main (int argc, char** argv) {
	try {
		return application{spam::config::from_args(argc, argv)}.exec();
	} catch (spam::config_exception const& e) {
		fmt::print(stderr, "{}", e.what());
		return 1;
	}
}
