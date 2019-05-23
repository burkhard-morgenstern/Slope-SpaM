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
