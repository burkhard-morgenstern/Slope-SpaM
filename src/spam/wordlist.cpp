#include "wordlist.hpp"

#include <algorithm>

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace fs = std::filesystem;
namespace rv = ranges::view;

namespace spam {

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
    words.erase(std::remove_if(words.begin(), words.end(),
        [](auto const& word) {
            return word == std::numeric_limits<word_t>::max();
        }),
        words.end());
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

} // namespace spam
