#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include <fmt/format.h>

#include "config.hpp"
#include "spam.hpp"

namespace fs = std::filesystem;

class application {
	spam::config config;
	std::shared_ptr<ThreadPool> threadpool;
	spam::pattern pattern;

public:
	application(spam::config config)
		: config(std::move(config)),
		threadpool(std::make_shared<ThreadPool>(
			std::thread::hardware_concurrency())),
		pattern(spam::pattern{this->config.pattern})
	{}

	auto exec()
		-> int
	{
		if (config.in.size() == 0) {
			fmt::print(stderr, "No input files!\n");
			return 1;
		}

		auto const pattern = spam::pattern{config.pattern};
		if (pattern.weight() > spam::wordlist::max_wordsize()) {
			fmt::print(stderr,
				"Unsupported pattern weight of {}! "
				"The supported maximum weight is {}!\n",
				pattern.weight(), spam::wordlist::max_wordsize());
			return 2;
		}

		auto threadpool = std::make_shared<ThreadPool>(
			std::thread::hardware_concurrency());
		if (config.out != "" && config.in.size() == 1) {
			auto outfile = std::ofstream{config.out};
			process(config.out, outfile);
			outfile.close();
		} else {
			for (auto& input_path : config.in) {
				auto output_path = !fs::is_directory(input_path)
					? fs::path{input_path}.replace_extension(".dmat")
					: fs::path{input_path.string() + ".dir.dmat"};
				auto file = std::ofstream{output_path};
				process(input_path, file);
				file.close();
			}
		}

		return 0;
	}

private:
	auto process(fs::path const& path, std::ostream& os)
		-> void
	{
		fmt::print("Processing {}...\n", path.string());
		fflush(stdout);
		auto sequences = fs::is_directory(path)
			? *spam::sequence::from_directory(path)
			: *spam::sequence::from_multi_fasta_file(path);
		if (sequences.size() == 0) {
			fmt::print(stderr, "Empty input path \"{}\".\n", path.string());
		} else {
			auto const matrix = spam::distance_matrix(
				std::move(sequences), pattern, threadpool);
			os << matrix;
		}
	}
};

int main (int argc, char** argv) {
	return application{spam::config::from_args(argc, argv)}.exec();
}
