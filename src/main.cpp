#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include "config.hpp"
#include "spam.hpp"

int main (int argc, char** argv) {
	auto config = spam::config::from_args(argc, argv);
	std::cout << "Loading sequences..." << std::flush;
	auto sequences = *spam::sequence::from_multi_fasta_file(config.in);
	std::cout << "\r" << sequences.size() << " sequences loaded!\n";
	auto pattern = spam::pattern{config.pattern};
	auto matrix = spam::distance_matrix(std::move(sequences), pattern);

	std::ofstream outfile(config.out);
	outfile << matrix;
	outfile.close();
}
