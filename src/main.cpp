#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include <fmt/format.h>
#include <range/v3/view.hpp>

#include "config.hpp"
#include "spam.hpp"

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
		for (auto const& [input_path, output_path] : input_output_mapping()) {
			auto file = std::ofstream{output_path};
			fmt::print("Processing {} -> {}{}",
				input_path.string(), output_path.string(), config.in.size() < 100 ? "\n" : "\r");
			process(input_path, file);
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

	auto process(fs::path const& path, std::ostream& os)
		-> void
	{
		fmt::print("Processing {}...\n", path.string());
		fflush(stdout);
		auto sequences = load_sequences(path);
		if (sequences.size() == 0) {
			fmt::print(stderr, "Empty input path \"{}\".\n", path.string());
		} else {
			auto wordlengths = config.wordlengths;
			if (wordlengths.empty()) {
				auto k = (spam::wordlist::kmin(sequences) +
					spam::wordlist::kmax(sequences)) / 2;
				fmt::print("No k provided! Using k = {}!\n", k);
				wordlengths.push_back(k);
			}
			auto const matrix = spam::distance_matrix(
				std::move(sequences),
				config.pattern,
				wordlengths,
				threadpool);
			os << matrix;
		}
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
