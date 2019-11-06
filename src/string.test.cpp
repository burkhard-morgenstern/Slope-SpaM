#include "catch2/catch.hpp"

#include "string.hpp"

TEST_CASE("split string") {
    SECTION("empty string") {
        REQUIRE(split("", ' ') == std::vector<std::string>{});
    }
    SECTION("simple") {
        REQUIRE(split("a,,b,", ',') ==
            std::vector<std::string>{"a", "", "b", ""});
    }
}
