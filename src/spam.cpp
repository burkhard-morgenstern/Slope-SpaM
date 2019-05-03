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
    spam::pattern const& pattern,
    size_t kmax)
    : _sequences(sequences),
    pattern(pattern),
    kmax(kmax),
    threadpool(std::thread::hardware_concurrency())
{
    calculate();
}

auto distance_matrix::sequences() const
    -> std::vector<spam::sequence> const&
{
    return _sequences;
}

auto distance_matrix::operator[](size_t idx) const
    -> std::vector<double> const&
{
    return matrix[idx];
}

void distance_matrix::calculate()
{
    initialize_matrix();
    create_wordlists_par();
    create_viewlists();
    calculate_matrix_par();
}

void distance_matrix::initialize_matrix()
{
    matrix.resize(sequences().size());
    for (auto& row : matrix) {
        row = std::vector<double>(sequences().size(), 0);
    }
}

auto create_spaced_words(
    spam::pattern const& pattern,
    std::string const& sequence)
    -> std::vector<std::string>
{
	auto result = std::set<std::string>{};
	for(auto substring = sequence.cbegin(); substring != sequence.cend() - pattern.bits.size() + 1; ++substring)
	{
		auto word = std::string{};
		word.reserve(pattern.size());
		for (auto idx : pattern) {
			word += substring[idx];
		}
		result.insert(std::move(word));
	}
	return {result.begin(), result.end()};
}

void distance_matrix::create_wordlists()
{
    wordlists = std::vector<std::vector<std::string>>{};
    wordlists.reserve(sequences().size());
    std::transform(sequences().begin(), sequences().end(),
        std::back_inserter(wordlists),
        [&](auto&& seq) { return create_spaced_words(pattern, seq.bases); });
}

void distance_matrix::create_wordlists_par() {
    wordlists = std::vector<std::vector<std::string>>{};
    wordlists.reserve(sequences().size());
    auto results = std::vector<std::future<std::vector<std::string>>>{};
    for (auto const& seq : sequences()) {
        results.push_back(threadpool.enqueue([&]() {
            return create_spaced_words(pattern, seq.bases);
        }));
    }
    auto n = results.size();
    auto i = 0;
    for (auto&& future : results) {
        std::cout << "\rCreating wordlists... " << i << "/" << n << std::flush;
        wordlists.push_back(future.get());
        i++;
    }
    std::cout << "\r";
}

void distance_matrix::create_viewlists()
{
    viewlists = std::vector<std::vector<std::pair<size_t, std::string_view>>>{};
    viewlists.reserve(wordlists.size());
    for (size_t i = 0; i < wordlists.size(); ++i) {
        auto viewlist = std::vector<std::pair<size_t, std::string_view>>{};
        viewlist.reserve(wordlists[i].size());
        std::transform(wordlists[i].begin(), wordlists[i].end(),
            std::back_inserter(viewlist),
            [&](auto const& word) {
                return std::make_pair(i, std::string_view{word});
            });
        viewlists.push_back(viewlist);
    }
}

void distance_matrix::calculate_matrix()
{
    for (size_t i = 0; i < sequences().size() - 1; i++) {
        for (size_t j = i + 1; j < sequences().size(); j++) {
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
    for (size_t i = 0; i < sequences().size() - 1; i++) {
        for (size_t j = i + 1; j < sequences().size(); j++) {
            results.push_back(
                threadpool.enqueue([&](size_t i, size_t j) {
                    return calculate_element(i, j);
                }, i, j)
            );
        }
    }

    size_t idx = 0;
    for (size_t i = 0; i < sequences().size() - 1; i++) {
        for (size_t j = i + 1; j < sequences().size(); j++) {
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
    int Lmax = std::max(sequences()[i].bases.size(), sequences()[j].bases.size());
    int kmin = log(2 * Lmax) / 0.87 + 1;
    auto wordlist = merge_viewlists(viewlists[i], viewlists[j]);
    auto const matches = calculate_matches(wordlist, kmin, kmax);
    return calculate_distance(matches, sequences()[i].bases.size(), sequences()[j].bases.size(), kmin);
}

auto distance_matrix::merge_viewlists(
    std::vector<std::pair<size_t, std::string_view>> const& lhs,
    std::vector<std::pair<size_t, std::string_view>> const& rhs)
    -> std::vector<std::pair<size_t, std::string_view>>
{
    auto wordlist = std::vector<std::pair<size_t, std::string_view>>{};
    wordlist.reserve(lhs.size() + rhs.size());
    std::merge(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(),
        std::back_inserter(wordlist),
        [](auto const& a, auto const& b) { return a.second < b.second; });
    return wordlist;
}

auto distance_matrix::calculate_matches(
    std::vector<std::pair<size_t, std::string_view>> wordlist,
    size_t kmin,
    size_t kmax)
    -> std::vector<size_t>
{
    auto matches = std::vector<size_t>{};
    matches.reserve(kmax - kmin);
    auto const i = wordlist[0].first;
    for (int k = kmin; k < kmax; k++)
    {
        size_t count = 0;
        auto it = wordlist.cbegin();
        while (it != wordlist.cend()) {
            auto const next = std::find_if_not(it, wordlist.cend(),
                [k, it](auto const& p) {
                    return std::equal(p.second.cbegin(),
                        p.second.cbegin() + k + 1, it->second.cbegin());
                });
            auto seq_i_count = std::count_if(it, next,
                [i](auto const& p) { return p.first == i; });
            auto seq_j_count = std::distance(it, next) - seq_i_count;
            count += seq_i_count * seq_j_count;
            it = next;
        }
        matches.push_back(count);
    }
    return matches;
}

auto distance_matrix::calculate_distance(std::vector<size_t> const& matches,
    size_t length1,
    size_t length2,
    size_t kmin)
    -> std::pair<double, double>
{
    double q = 0.25;
    std::vector<long double> x_values;
    std::vector<long double> y_values;
    x_values.reserve(matches.size());
    y_values.reserve(matches.size());
    for (int w = 0; w < matches.size(); w++){
        if (matches[w] != 0){
            long double x = pow(q, kmin + w + 1) * length1 * length2; // TODO:again, why + 1
            long double y = log(matches[w] - x);
            if (std::isnan(y) == false) {
                x_values.push_back(kmin + w + 1);
                y_values.push_back(y);
            }
        }
    }
    auto m = *slope(
        x_values.begin(), x_values.end(),
        y_values.begin(), y_values.end());
    double p = exp(m);
    double d = -(3.0 / 4.0) * log(1 - (4.0 / 3.0) * (1 - p));
    return {p, d};
}

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix) {
	auto const& sequences = matrix.sequences();
	os << sequences.size() << std::endl;
	for (int i = 0; i < sequences.size(); i++){
		os << sequences[i].name << '\t';
		for (int k = 0; k < sequences.size(); k++){
			os << matrix[i][k] << '\t';
		}
		os << std::endl;
	}
	return os;
}

} // namespace spam
