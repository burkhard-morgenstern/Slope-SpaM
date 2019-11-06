#pragma once

#include <filesystem>
#include <regex>

auto resolve_wildcards(std::string const& wildcard_path)
	-> std::vector<std::filesystem::path>;
