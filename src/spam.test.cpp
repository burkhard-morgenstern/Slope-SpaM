#include "catch2/catch.hpp"

#include "spam.hpp"

TEST_CASE("pattern") {
    SECTION("simple pattern") {
        auto pattern = spam::pattern{"1111"};
        REQUIRE(pattern.weight() == 4);
        auto expected = size_t{0};
        for (auto v : pattern) {
            REQUIRE(v == expected);
            expected++;
        }
    }
    SECTION("complex pattern") {
        auto pattern = spam::pattern{"1010101"};
        REQUIRE(pattern.weight() == 4);
        auto expected = size_t{0};
        for (auto v : pattern) {
            REQUIRE(v == expected);
            expected += 2;
        }
    }
    SECTION("reduce"){
        SECTION("invalid k") {
            auto pattern = spam::pattern{"1111"};
            REQUIRE_THROWS_AS(pattern.reduce(5), std::invalid_argument);
        }
        SECTION("no-op") {
            auto pattern = spam::pattern{"1111"};
            REQUIRE(pattern.reduce(4) == pattern);
        }
        SECTION("valid k, simple") {
            auto pattern = spam::pattern{"1111"}.reduce(3);
            REQUIRE(pattern == spam::pattern{"111"});
        }
        SECTION("valid k, complex") {
            auto pattern = spam::pattern{"10111"}.reduce(3);
            REQUIRE(pattern == spam::pattern{"1011"});
            REQUIRE(pattern.size() == 4);
            REQUIRE(pattern.weight() == 3);
        }
        SECTION("trailing zeros") {
            auto a = spam::pattern{"1011"};
            auto b = spam::pattern{"10110"};
            REQUIRE(a == b);
            REQUIRE(a.size() == b.size());
            REQUIRE(a.weight() == b.weight());
        }
    }
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
        auto expected = std::vector<size_t>{0x1B, 0x6C, 0xB1, 0xC6};
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
    SECTION("complex pattern") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", "AACCGGTTAA"}};
        auto pattern = spam::pattern{"1010101"};
        auto expected = std::vector<size_t>{0x1B, 0x1B, 0x6C, 0x6C};
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
    SECTION("sorting") {
        auto sequence = spam::sequence{spam::assembled_sequence{"", "TGCATGCATGC"}};
        auto pattern = spam::pattern{"1111"};
        auto expected = std::vector<size_t>{
            0x39, 0x39, 0x4E, 0x4E, 0x93, 0x93, 0xE4, 0xE4};
        auto wordlist = spam::wordlist{sequence, pattern};
        REQUIRE(std::equal(wordlist.begin(), wordlist.end(), expected.begin()));
    }
}
