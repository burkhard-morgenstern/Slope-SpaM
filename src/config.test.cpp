#include "catch2/catch.hpp"

#include "config.hpp"
#include "spam/wordlist.hpp"

namespace spam {

auto parse_pattern(std::string const& patternflag)
    -> spam::pattern;

} // namespace spam

TEST_CASE("patternflag parser") {
    SECTION("empty patternflag") {
        REQUIRE(spam::parse_pattern("").weight() == 0);
    }
    SECTION("simple pattern") {
        auto pattern = spam::parse_pattern("110011");
        REQUIRE(pattern.size() == 6);
        REQUIRE(pattern.weight() == 4);
        auto indices = std::vector<size_t>{0, 1, 4, 5};
        REQUIRE(std::equal(
            pattern.begin(), pattern.end(),
            indices.begin(), indices.end()));
    }
    SECTION("random pattern") {
        SECTION("1.0 probability") {
            auto pattern = spam::parse_pattern("1.0");
            REQUIRE(pattern.size() == spam::wordlist::max_wordsize());
            REQUIRE(pattern.weight() == spam::wordlist::max_wordsize());
        }
        SECTION("0.0 probablility") {
            REQUIRE_THROWS_AS(spam::parse_pattern("0.0"), spam::config_exception);
        }
        SECTION("0.5 probablility") {
            REQUIRE(spam::parse_pattern("0.5").weight()
                == spam::wordlist::max_wordsize());
        }
    }
}
