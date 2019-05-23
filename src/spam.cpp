#include "spam.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <thread>

#include "math.hpp"

namespace spam {

pattern::pattern(std::string _bits)
    : bits(std::move(_bits))
{
    auto end = static_cast<ptrdiff_t>(bits.size());
    for (; end >= 0 && bits[end] != '1'; --end);
    bits.resize(end + 1);
    for (size_t i = 0; i < bits.size(); ++i) {
        if (bits[i] == '1') {
            indices.push_back(i);
        }
    }
}

auto pattern::reduce(size_t k) const
    -> pattern
{
    if (weight() < k) {
        throw std::invalid_argument{"weight < k"};
    }

    if (weight() == k) {
        return *this;
    }

    return pattern{{bits.begin(), bits.begin() + indices[k]}};
}

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
    std::string bases;
    std::getline(is, bases, '>');
    bases.erase(std::remove_if(bases.begin(), bases.end(),
        [](const char c) { return isspace(c); }), bases.end());
    seq.bases = std::move(bases);
    return is;
}

auto encode_sequence(spam::sequence const& sequence)
    -> std::vector<int>
{
    auto result = std::vector<int>{};
    result.reserve(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
        std::back_inserter(result),
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
    return result;
}

template<class EncodedIt>
auto create_word(EncodedIt it, spam::pattern const& pattern)
    -> word_t
{
    auto word = word_t{0};
    auto successful = true;
    for (auto pos : pattern) {
        if (it[pos] == -1) {
            successful = false;
            break;
        }
        word <<= 2;
        word += it[pos];
    }
    return successful
        ? word
        : std::numeric_limits<word_t>::max();
}

wordlist::wordlist(
    spam::sequence const& sequence,
    spam::pattern pattern)
    : pattern(pattern)
{
    if (sequence.size() < pattern.size()) {
        return;
    }
    auto encoded = encode_sequence(sequence);
    words.reserve(sequence.size() - pattern.size() + 1);
    for (auto it = encoded.begin();
        it != encoded.end() - pattern.size() + 1;
        ++it)
    {
        words.push_back(create_word(it, pattern));
    }
    std::sort(words.begin(), words.end());
}

distance_matrix::distance_matrix(
    std::vector<spam::sequence>&& sequences,
    spam::pattern pattern)
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
    create_wordlists();
    calculate_matrix();
}

void distance_matrix::initialize_matrix()
{
    matrix.resize(sequences.size());
    for (auto& row : matrix) {
        row = std::vector<double>(sequences.size(), 0);
    }
}

void distance_matrix::create_wordlists()
{
    auto Lmax = std::max_element(
        sequences.begin(), sequences.end(),
        [](auto&& a, auto&& b) {
            return a.size() < b.size();
        })->size();
    size_t kmin = std::ceil(std::log(Lmax) / 0.87);
    size_t kmax = std::floor(std::log(Lmax) / 0.63);
    auto k = (kmin + kmax) / 2;
    pattern = pattern.reduce(k);
    wordlists.reserve(sequences.size());
    auto results = std::vector<std::future<spam::wordlist>>{};
    for (auto i = size_t{0}; i < sequences.size(); ++i) {
        results.push_back(threadpool.enqueue([&, i = i]() {
            return spam::wordlist(sequences[i], pattern);
        }));
    }
    for (auto&& future : results) {
        wordlists.push_back(future.get());
    }
}

void distance_matrix::calculate_matrix()
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
            matrix[i][j] = distance;
            matrix[j][i] = distance;
        }
    }
}

auto calculate_matches(
    spam::wordlist const& wordlist1,
    spam::wordlist const& wordlist2)
    -> size_t
{
    size_t count = 0;
    auto next_fn = [](auto it, auto&& wordlist) {
        return std::find_if_not(it, wordlist.end(),
            [it](auto const& p) {
                return p == *it;
            });
    };
    auto it1 = wordlist1.begin();
    auto next1 = next_fn(it1, wordlist1);
    auto it2 = wordlist2.begin();
    auto next2 = next_fn(it2, wordlist2);
    while (it1 != wordlist1.end() && it2 != wordlist2.end()) {
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

auto calculate_distance(
    size_t matches,
    spam::pattern const& pattern,
    spam::sequence const& seq1,
    spam::sequence const& seq2)
    -> std::pair<double, double>
{
    auto length1 = seq1.size() - pattern.weight() + 1;
    auto length2 = seq2.size() - pattern.weight();
    auto q = 0.25;
    long double e = pow(q, pattern.weight()) * length1 * length2;
    long double y = log(matches - e) - log(length1);
    long double x = pattern.weight();
    auto m = y / x;
    auto p = exp(m);
    auto d = -(3.0 / 4.0) * log(1 - (4.0 / 3.0) * (1 - p));
    return {p, d};
}

auto distance_matrix::calculate_element(size_t i, size_t j) const
    -> std::pair<double, double>
{
    return calculate_distance(
        calculate_matches(wordlists[i], wordlists[j]),
        pattern, sequences[i], sequences[j]);
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
