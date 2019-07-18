#include "spam.hpp"

#include <algorithm>
#include <fstream>
#include <set>
#include <thread>

#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include "math.hpp"

namespace fs = std::filesystem;
namespace rv = ranges::view;

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

auto load_fasta_file(std::filesystem::path const& filename)
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

auto wordlist::max_wordsize()
    -> size_t
{
    return 4 * sizeof(word_t) - 1;
}

auto wordlist::kmin(sequence const& seq)
    -> size_t
{
    return std::ceil(std::log(seq.size()) / 0.87);
}

auto wordlist::kmin(std::vector<sequence> const& seqs)
    -> size_t
{
    return ranges::max(
        seqs | rv::transform([](auto&& seq) { return kmin(seq); }));
}

auto wordlist::kmax(sequence const& seq)
    -> size_t
{
    return std::floor(std::log(seq.size()) / 0.63);
}

auto wordlist::kmax(std::vector<sequence> const& seqs)
    -> size_t
{
    return ranges::min(
        seqs | rv::transform([](auto&& seq) { return kmax(seq); }));
}

auto encode_sequence(std::string const& sequence)
    -> std::vector<int8_t>
{
    auto result = std::vector<int8_t>{};
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
    word <<= 8 * sizeof(word_t) - 2 * pattern.weight();
    return successful
        ? word
        : std::numeric_limits<word_t>::max();
}

auto create_words(std::string const& sequence, spam::pattern const& pattern)
    -> std::vector<word_t>
{
    if (sequence.size() < pattern.size()) {
        return {};
    }
    auto encoded = encode_sequence(sequence);
    auto words = std::vector<word_t>{};
    words.reserve(sequence.size() - pattern.size() + 1);
    for (auto it = encoded.begin();
        it != encoded.end() - pattern.size() + 1;
        ++it)
    {
        words.push_back(create_word(it, pattern));
    }
    return words;
}

wordlist::wordlist(
    spam::sequence const& sequence,
    spam::pattern pattern)
    : pattern(pattern)
{
    if (std::holds_alternative<assembled_sequence>(sequence)) {
        words = create_words(
            std::get<assembled_sequence>(sequence).nucleotides, pattern);
    } else {
        auto& reads = std::get<unassembled_sequence>(sequence).reads;
        for (auto& read : reads) {
            auto tmp = create_words(read, pattern);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(words));
        }
    }
    std::sort(words.begin(), words.end());
}

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
        wordlists.emplace_back(sequences[i], pattern);
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

auto reduce_wordlist(spam::wordlist const& wordlist, size_t k)
{
    auto reduction_pattern = std::numeric_limits<word_t>::max();
    reduction_pattern <<= 8 * sizeof(word_t) - 2 * k;
    return wordlist
        | rv::transform(
            [&wordlist, reduction_pattern = reduction_pattern](auto word) {
                return word & reduction_pattern;
            });
}

auto calculate_matches(
    spam::wordlist const& _wordlist1,
    spam::wordlist const& _wordlist2,
    size_t k)
    -> size_t
{
    auto wordlist1 = reduce_wordlist(_wordlist1, k);
    auto wordlist2 = reduce_wordlist(_wordlist2, k);

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
            count += std::min(
                std::distance(it1, next1),
                std::distance(it2, next2));
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
    size_t kmax,
    spam::sequence const& seq1,
    spam::sequence const& seq2)
    -> std::pair<double, double>
{
    auto q = background_match_probability(seq1, seq2);
    auto lh = std::min(seq1.size() + 1, seq2.size());
    auto values = matches
        | rv::transform(
            [&](auto&& p) {
                auto [k, matches] = p;
                return std::pair<double, double>(
                    k, log(matches) - log(lh - kmax));
            })
        | ranges::to<std::vector<std::pair<double, double>>>();
    if (values.size() == 1) {
        values.emplace_back(0, 0);
    }
    auto m = slope(values);
    auto p = exp(m);
    auto d = -(3.0 / 4.0) * log(1 - (4.0 / 3.0) * (1 - p));
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
    return calculate_distance(
        matches, wordlist::kmax(sequences),
        sequences[i], sequences[j]);
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
