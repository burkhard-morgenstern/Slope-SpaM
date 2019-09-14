#include "distance_matrix.hpp"

#include <thread>

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include "../math.hpp"

namespace fs = std::filesystem;
namespace rv = ranges::view;

namespace spam {

distance_matrix::distance_matrix(
    std::vector<spam::sequence> sequences,
    spam::pattern pattern,
    std::vector<size_t> wordlengths,
    std::optional<std::shared_ptr<ThreadPool>> threadpool /* = {} */)
    : sequences(std::move(sequences)),
    pattern(pattern),
    wordlengths(std::move(wordlengths)),
    threadpool(std::move(threadpool))
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
    auto max_wordlength = *ranges::max_element(wordlengths);
    pattern = pattern.reduce(max_wordlength);
    if (!threadpool) {
        create_wordlists_seq();
    } else {
        create_wordlists_par();
    }
}

void distance_matrix::create_wordlists_seq()
{
    wordlists.reserve(sequences.size());
    auto results = std::vector<std::future<spam::wordlist>>{};
    for (auto i = size_t{0}; i < sequences.size(); ++i) {
        wordlists.push_back(spam::wordlist(sequences[i], pattern));
    }
}

void distance_matrix::create_wordlists_par()
{
    wordlists.reserve(sequences.size());
    auto results = std::vector<std::future<spam::wordlist>>{};
    for (auto i = size_t{0}; i < sequences.size(); ++i) {
        results.push_back((*threadpool)->enqueue([&, i = i]() {
            return spam::wordlist(sequences[i], pattern);
        }));
    }
    for (auto&& future : results) {
        wordlists.push_back(future.get());
    }
}

void distance_matrix::calculate_matrix()
{
    if (!threadpool) {
        calculate_matrix_seq();
    } else {
        calculate_matrix_par();
    }
}

void distance_matrix::calculate_matrix_seq()
{
    for (size_t i = 0; i < sequences.size() - 1; i++) {
        for (size_t j = i + 1; j < sequences.size(); j++) {
            auto const [probability, distance] = calculate_element(i, j);
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
                (*threadpool)->enqueue([&](size_t i, size_t j) {
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
    spam::wordlist const& _wordlist1,
    spam::wordlist const& _wordlist2,
    size_t k)
    -> size_t
{
    auto wordlist1 = _wordlist1.reduce(k);
    auto wordlist2 = _wordlist2.reduce(k);

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
        auto const v1 = *it1;
        auto const v2 = *it2;
        if (v1 == v2) {
            count += std::distance(it1, next1) * std::distance(it2, next2);
        }
        if (v1 <= v2) {
            it1 = next1;
            next1 = next_fn(it1, wordlist1);
        }
        if (v1 >= v2) {
            it2 = next2;
            next2 = next_fn(it2, wordlist2);
        }
    }
    return count;
}

auto count_nucleotide(spam::sequence const& sequence, char nucleotide)
    -> size_t
{
    if (std::holds_alternative<assembled_sequence>(sequence)) {
        return ranges::count(std::get<assembled_sequence>(sequence).nucleotides, nucleotide);
    } else {
        return ranges::accumulate(
            std::get<unassembled_sequence>(sequence).reads
                | rv::transform(
                    [nucleotide](auto&& read) {
                        return ranges::count(read, nucleotide);
                    }),
            size_t{0});
    }
}

auto background_match_probability(
    spam::sequence const& seq1,
    spam::sequence const& seq2)
    -> double
{
    auto result = 0.0;
    for (auto c : {'A', 'C', 'G', 'T'}) {
        result +=
            1.0 * count_nucleotide(seq1, c) / seq1.size()
            * (1.0 * count_nucleotide(seq2, c) / seq2.size());
    }
    return result;
}

template<class Range>
auto calculate_distance(
    Range&& matches,
    spam::sequence const& seq1,
    spam::sequence const& seq2)
    -> std::pair<double, double>
{
    auto q = 0.25;
    auto values = matches
        | rv::transform(
            [&](auto&& p) {
                auto [k, matches] = p;
                long double e = pow(q, k) * seq1.adjusted_size(k) * seq2.adjusted_size(k);
                return std::pair<long double, long double>(k, log(matches - e));
            })
        | ranges::to<std::vector<std::pair<long double, long double>>>();
    auto m = slope(values);
    auto p = exp(m) / ((1.0 - seq1.error_rate()) * (1.0 - seq2.error_rate()));
    auto d = -(3.0 / 4.0) * log(1.0 - (4.0 / 3.0) * (1.0 - p));
    return {p, d};
}

auto distance_matrix::calculate_element(size_t i, size_t j) const
    -> std::pair<double, double>
{
    auto matches = wordlengths
        | rv::transform(
            [&](auto k) {
                return std::make_pair(
                    k, calculate_matches(wordlists[i], wordlists[j], k));
            });
    auto kmax = *ranges::max_element(wordlengths);
    return calculate_distance(matches, sequences[i], sequences[j]);
}

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix)
{
	os << matrix.size() << std::endl;
	for (int i = 0; i < matrix.size(); i++) {
        auto const& [sequence, distances] = matrix.column(i);
		os << sequence.name() << '\t';
		for (int j = 0; j < matrix.size(); j++){
			os << distances[j] << '\t';
		}
		os << "\n";
	}
	return os;
}

} // namespace spam
