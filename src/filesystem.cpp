#include "filesystem.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

namespace fs = std::filesystem;
namespace rv = ranges::view;

auto split(std::string const& s, char seperator)
	-> std::vector<std::string>
{
	auto output = std::vector<std::string>{};
	auto prev_pos = size_t{0};
	auto pos = size_t{0};
	while((pos = s.find(seperator, pos)) != std::string::npos) {
		auto substring = s.substr(prev_pos, pos - prev_pos);
		output.push_back(substring);
		prev_pos = ++pos;
	}
	output.push_back(s.substr(prev_pos, pos - prev_pos));
	return output;
}

auto regex_escape(std::string const& regex_string)
    -> std::string
{
    auto result = std::string{"^"};
    for (auto c : regex_string) {
        switch (c) {
        case '.':
            result += "\\.";
            break;
        case '*':
            result += ".*";
            break;
        default:
            result += c;
            break;
        }
    }
    result += "$";
    return result;
}

auto components(std::filesystem::path const& path)
{
    auto path_str = std::string{path};
    return path_str
		| ranges::view::split('/')
		| ranges::view::transform(
			[](auto&& rng) {
            	return std::string{
                    std::string_view(&*rng.begin(), ranges::distance(rng))};
			})
		| ranges::view::remove_if(
			[](auto&& str) { return str == ""; })
        | ranges::to<std::vector<std::string>>();
}

auto children(std::filesystem::path const& directory)
{
    return std::filesystem::directory_iterator(directory)
        | ranges::view::transform(
            [](auto&& entry) {
                return entry.path();
            });
}

auto matching_children(
    std::filesystem::path const& directory,
    std::regex const& component_regex)
{
    return children(directory)
        | ranges::view::remove_if(
            [&component_regex](auto&& path) { 
                return !std::regex_match(
                    path.string(), component_regex);
            });
}

template<class Range>
auto all_matching_children(
    Range&& directories,
    std::regex const& component_regex)
    -> std::vector<fs::path>
{
    return directories
        | rv::transform(
            [&component_regex](auto&& path) {
                return matching_children(path, component_regex);
            })
        | rv::join
        | ranges::to<std::vector<fs::path>>();
}

template<class Range>
auto all_matching_children(
    Range&& directories,
    std::string const& component)
    -> std::vector<fs::path>
{
    return directories
        | rv::transform(
            [&component](auto&& path) {
                return path / component;
            })
        | rv::remove_if(
            [](auto&& path) {
                return !fs::exists(path);;
            })
        | ranges::to<std::vector<fs::path>>();
}

auto subdirectories(std::filesystem::path const& directory)
{
    return children(directory)
        | ranges::view::remove_if(
            [](auto&& path) {
                return !std::filesystem::is_directory(path);
            });
}

template<class Range>
auto all_subdirectories(Range&& paths)
    -> std::vector<fs::path>
{
    return paths
        | rv::transform(
            [](auto&& path) {
                return subdirectories(path);
            })
        | rv::join
        | ranges::to<std::vector<fs::path>>();
}

auto resolve_wildcards(std::string const& wildcard_path)
	-> std::vector<fs::path>
{
	if (ranges::count(wildcard_path, '*') == 0) {
		return {wildcard_path};
	}
	auto const _components = components(wildcard_path);
	if (ranges::empty(_components)) {
		return {};
	}
	auto paths = std::vector<fs::path>{"."};
	for (auto const component : _components) {
        auto existing_paths = paths
			| rv::remove_if([](auto&& path) { return !fs::exists(path); });
        auto const count = ranges::count(component, '*');
        if (count == 2 && component == "**") {
            paths = all_subdirectories(existing_paths);
        } else {
            if (count > 0) {
		        auto const regex = std::regex{regex_escape(component)};
                paths = all_matching_children(existing_paths, regex);
            } else {
                paths = all_matching_children(existing_paths, component);
            }
        }
	}
	return paths;
}
