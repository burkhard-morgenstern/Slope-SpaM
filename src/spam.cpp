#include "spam.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <thread>

#include "math.hpp"

namespace spam {

auto sequence::from_multi_fasta_file(std::string const& filename)
    -> std::optional<std::vector<sequence>>
{
    auto file = std::ifstream(filename);
    if (!file.is_open()) {
        return {};
    }
	auto result = std::vector<spam::sequence>{};
	auto dummy = std::string{};
    std::getline(file, dummy, '>');
    while (!file.eof()) {
		auto seq = spam::sequence{};
		file >> seq;
		result.push_back(std::move(seq));
	}
	return {result};
}

auto operator>>(std::istream& is, sequence& seq)
    -> std::istream&
{
    std::getline(is, seq.name);
    std::getline(is, seq.bases, '>');
    seq.bases.erase(std::remove_if(seq.bases.begin(), seq.bases.end(),
        [](const char c) { return isspace(c); }), seq.bases.end());
    return is;
}

distance_matrix::distance_matrix(
    std::vector<spam::sequence> const& sequences,
    spam::pattern const& pattern)
    : distance_matrix(std::vector<spam::sequence>{sequences}, pattern)
{}

distance_matrix::distance_matrix(
    std::vector<spam::sequence>&& sequences,
    spam::pattern const& pattern)
    : sequences(std::move(sequences)),
    pattern(pattern),
    threadpool(std::thread::hardware_concurrency())
{
    calculate();
}

auto distance_matrix::size() const
    -> size_t
{
    return sequences.size();
}

auto distance_matrix::column(size_t i) const
    -> std::pair<spam::sequence const&, std::vector<double> const&>
{
    return {sequences[i], matrix[i]};
}

void distance_matrix::calculate()
{
    initialize_matrix();
    create_wordlists_par();
    calculate_matrix_par();
}

void distance_matrix::initialize_matrix()
{
    matrix.resize(sequences.size());
    for (auto& row : matrix) {
        row = std::vector<double>(sequences.size(), 0);
    }
}

auto create_spaced_words(
    spam::pattern const& pattern,
    std::string const& sequence,
    size_t k)
    -> std::vector<word_t>
{
    auto seq_values = std::vector<int>{};
    seq_values.reserve(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
        std::back_inserter(seq_values),
        [](auto&& c) {
            switch(c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            }
            return -1;
        });
    seq_values.erase(std::remove_if(
        seq_values.begin(), seq_values.end(),
        [](auto&& v) { return v == -1; }), seq_values.end());

	auto result = std::vector<word_t>{};
    result.reserve(sequence.size() - k + 1);
	for(auto substring = seq_values.begin();
        substring != seq_values.end() - k + 1;
        ++substring)
	{
        auto word = word_t{0};
        for (auto i = size_t{0}; i < k; ++i) {
            word <<= 2;
            word += substring[pattern[i]];
        }
        result.push_back(word);
	}
    std::sort(result.begin(), result.end());
	return result;
}

void distance_matrix::create_wordlists_par()
{
    auto Lmax = std::max_element(
        sequences.begin(), sequences.end(),
        [](auto&& a, auto&& b) {
            return a.bases.size() < b.bases.size();
        })->bases.size();
    size_t kmin = std::ceil(std::log(Lmax) / 0.87);
    size_t kmax = std::floor(std::log(Lmax) / 0.63);
    k = (kmin + kmax) / 2;
    wordlists.reserve(sequences.size());
    auto results = std::vector<std::future<std::vector<word_t>>>{};
    for (auto i = size_t{0}; i < sequences.size(); ++i) {
        results.push_back(threadpool.enqueue([&, i = i]() {
            return create_spaced_words(pattern, sequences[i].bases, k);
        }));
    }
    auto n = results.size();
    auto i = 0;
    for (auto&& future : results) {
        std::cout << "\rCreating wordlists... " << i << "/" << n << std::flush;
        wordlists.push_back(future.get());
        ++i;
    }
    std::cout << "\rFinished...";
}

void distance_matrix::calculate_matrix()
{
    for (size_t i = 0; i < sequences.size() - 1; i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            auto const [probability, distance] = calculate_element(i, j);
            std::cout << "match probability p : " << probability << " Jukes-Cantor distance d : " << distance << '\n';
            matrix[i][j] = distance;
            matrix[j][i] = distance;
        }
    }
}

void distance_matrix::calculate_matrix_par()
{
    auto results = std::vector<std::future<std::pair<double, double>>>{};
    for (size_t i = 0; i < sequences.size() - 1; i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            results.push_back(
                threadpool.enqueue([&](size_t i, size_t j) {
                    return calculate_element(i, j);
                }, i, j)
            );
        }
    }

    size_t idx = 0;
    for (size_t i = 0; i < sequences.size() - 1; i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            auto const [probability, distance] = results[idx++].get();
            std::cout << "match probability p : " << probability << " Jukes-Cantor distance d : " << distance << '\n';
            matrix[i][j] = distance;
            matrix[j][i] = distance;
        }
    }
}

auto distance_matrix::calculate_element(size_t i, size_t j) const
    -> std::pair<double, double>
{
    auto const matches = calculate_matches(wordlists[i], wordlists[j]);
    return calculate_distance(
        matches,
        sequences[i].bases.size() - k + 1,
        sequences[j].bases.size() - k,
        k);
}

auto distance_matrix::calculate_matches(
		std::vector<word_t> const& wordlist1,
		std::vector<word_t> const& wordlist2)
    -> size_t
{
    size_t count = 0;
    auto next_fn = [](auto it, auto&& wordlist) {
        return std::find_if_not(it, wordlist.cend(),
            [it](auto const& p) {
                return p == *it;
            });
    };
    auto it1 = wordlist1.begin();
    auto next1 = next_fn(it1, wordlist1);
    auto it2 = wordlist2.begin();
    auto next2 = next_fn(it2, wordlist2);
    while (it1 != wordlist1.end() && it2 != wordlist2.cend()) {
        if (*it1 == *it2) {
            count += std::distance(it1, next1) * std::distance(it2, next2);
            it1 = next1;
            next1 = next_fn(it1, wordlist1);
            it2 = next2;
            next2 = next_fn(it2, wordlist2);
        } else {
            if (*it1 < *it2) {
                it1 = next1;
                next1 = next_fn(it1, wordlist1);
            } else {
                it2 = next2;
                next2 = next_fn(it2, wordlist2);
            }
        }
    }
    return count;
}

auto distance_matrix::calculate_distance(
    size_t matches,
    size_t length1,
    size_t length2,
    size_t k)
    -> std::pair<double, double>
{
    double q = 0.25;
    long double x = pow(q, k) * length1 * length2;
    long double y = log(matches - x) - log(length1);
    x = k;
    auto m = y / x;
    auto p = exp(m);
    auto d = -(3.0 / 4.0) * log(1 - (4.0 / 3.0) * (1 - p));
    return {p, d};
}

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix)
{
	os << matrix.size() << std::endl;
	for (int i = 0; i < matrix.size(); i++) {
        auto const& [sequence, distances] = matrix.column(i);
		os << sequence.name << '\t';
		for (int j = 0; j < matrix.size(); j++){
			os << distances[j] << '\t';
		}
		os << "\n";
	}
	return os;
}

} // namespace spam
