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

std::vector<spam::sequence> read_sequences(std::string const& file_name)
{
    std::ifstream infile(file_name);
    if (infile.is_open() == false)
    {
        std::cerr << "File " + file_name + " could not be opened!";
		exit(-1);
    }
    std::string line;
    std::string header;
    std::getline(infile, line, '>');
	auto result = std::vector<spam::sequence>{};
    while (!infile.eof())
    {
        std::getline(infile, header);
        std::getline(infile, line, '>');
		line.erase(std::remove_if(line.begin(), line.end(),
			[](const char c) { return isspace(c); }), line.end());
		result.emplace_back(spam::sequence{header, line});
	}
	std::cout << result.size() << " sequences have been read." << std::endl;
	return result;
}

int check_pattern(spam::pattern const& pattern, std::vector<spam::sequence> const& sequences)
{
	int lenPat = pattern.bits.size();
	int kmax = pattern.indices.size();
	int Lm  = std::max_element(sequences.begin(), sequences.end(),
		[](auto const& a, auto const& b) { return a.bases.size() < b.bases.size(); })->bases.size();
	int km =  log (2 * Lm ) / 0.87 + 1;
	if (kmax-km < 10){
		// TODO: was does this even mean?
		std::cout << "Your pattern weight is too small for your sequences, there should be a least 10 value calculable. Please use a pattern with a higher pattern weight."<< '\n';
		exit(-1);
	}
	return kmax;
}

int main (int argc, char** argv){
	auto config = spam::config::from_args(argc, argv);
	auto sequences = read_sequences(config.in);
	auto pattern = spam::pattern{config.pattern};
	int kmax = check_pattern(pattern, sequences);

	auto matrix = distance_matrix(sequences, pattern, kmax);
	std::ofstream file(config.out);
	file << matrix;
	file.close();
}
