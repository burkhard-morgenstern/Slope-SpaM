#include "catch2/catch.hpp"

#include "wordlist.hpp"

template<class Range>
auto shifted_words(Range&& words, size_t pattern_weight)
    -> std::vector<spam::word_t>
{
    auto result = std::vector<spam::word_t>{};
    result.reserve(words.size());
    std::transform(words.begin(), words.end(),
        std::back_inserter(result),
        [pattern_weight = pattern_weight](auto&& word) {
            return word << (8 * sizeof(spam::word_t) - 2 * pattern_weight);
        });
    return result;
}

TEST_CASE("wordlist") {
    SECTION("insufficient sequence length") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", ""}};
        auto pattern = spam::pattern{"11111111"};
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(wordlist.begin() == wordlist.end());
    }
    SECTION("simple pattern") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", "ACGTACG"}};
        auto pattern = spam::pattern{"1111"};
        auto expected = shifted_words(std::vector<spam::word_t>{0x1B, 0x6C, 0xB1, 0xC6}, 4);
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
    SECTION("complex pattern") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", "AACCGGTTAA"}};
        auto pattern = spam::pattern{"1010101"};
        auto expected = shifted_words(std::vector<spam::word_t>{0x1B, 0x1B, 0x6C, 0x6C}, 4);
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
    SECTION("sorting") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", "TGCATGCATGC"}};
        auto pattern = spam::pattern{"1111"};
        auto expected = shifted_words(std::vector<spam::word_t>{
            0x39, 0x39, 0x4E, 0x4E, 0x93, 0x93, 0xE4, 0xE4}, 4);
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
}

TEST_CASE("wordlist matches") {
    auto sequence = spam::sequence{spam::assembled_sequence{"", "TGCATGCATCG"}};
    auto pattern = spam::pattern{"1111"};
    auto wordlist = spam::wordlist{sequence, pattern};
    REQUIRE(calculate_matches(wordlist, wordlist, 4) == 12);
    REQUIRE(calculate_matches(wordlist, wordlist, 3) == 14);
    REQUIRE(calculate_matches(wordlist, wordlist, 2) == 16);
}
