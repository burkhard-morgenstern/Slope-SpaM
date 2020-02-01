#include "filesystem.hpp"

#include "string.hpp"

namespace fs = std::filesystem;

auto regex_escape(std::string const& regex_string) -> std::string {
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

auto components(std::filesystem::path const& path) -> std::vector<std::string> {
  auto components = split(path, '/');
  components.erase(
      std::remove_if(components.begin(), components.end(),
                     [](auto& component) { return component == ""; }),
      components.end());
  return components;
}

auto children(std::filesystem::path const& directory)
    -> std::vector<std::filesystem::path> {
  auto result = std::vector<fs::path>{};
  for (auto& entry : std::filesystem::directory_iterator(directory)) {
    result.push_back(entry.path());
  }
  return result;
}

auto matching_children(std::filesystem::path const& directory,
                       std::regex const& component_regex)
    -> std::vector<std::filesystem::path> {
  auto result = children(directory);
  result.erase(std::remove_if(result.begin(), result.end(),
                              [&component_regex](auto& path) {
                                return !std::regex_match(path.string(),
                                                         component_regex);
                              }),
               result.end());
  return result;
}

template <class Range>
auto all_matching_children(Range&& directories,
                           std::regex const& component_regex)
    -> std::vector<fs::path> {
  auto result = std::vector<fs::path>{};
  for (auto& directory : directories) {
    auto matching = matching_children(directory, component_regex);
    std::copy(matching.begin(), matching.end(), std::back_inserter(result));
  }
  return result;
}

template <class Range>
auto all_matching_children(Range&& directories, std::string const& component)
    -> std::vector<fs::path> {
  auto result = std::vector<fs::path>{};
  result.reserve(directories.size());
  std::transform(
      directories.begin(), directories.end(), std::back_inserter(result),
      [&component](auto&& directory) { return directory / component; });
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](auto&& path) { return !fs::exists(path); }),
               result.end());
  return result;
}

auto subdirectories(std::filesystem::path const& directory) {
  auto result = children(directory);
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](auto&& path) {
                                return !std::filesystem::is_directory(path);
                              }),
               result.end());
  return result;
}

template <class Range>
auto all_subdirectories(Range&& paths) -> std::vector<fs::path> {
  auto result = std::vector<fs::path>{};
  for (auto& path : paths) {
    auto subdirs = subdirectories(path);
    std::copy(subdirs.begin(), subdirs.end(), std::back_inserter(result));
  }
  return result;
}

auto resolve_wildcards(std::string const& wildcard_path)
    -> std::vector<fs::path> {
  if (std::count(wildcard_path.begin(), wildcard_path.end(), '*') == 0) {
    return {wildcard_path};
  }

  auto const _components = components(wildcard_path);
  if (std::empty(_components)) {
    return {};
  }

  auto paths = std::vector<fs::path>{"."};
  for (auto const component : _components) {
    auto existing_paths = paths;
    existing_paths.erase(
        std::remove_if(existing_paths.begin(), existing_paths.end(),
                       [](auto&& path) { return !fs::exists(path); }),
        existing_paths.end());
    auto const count = std::count(component.begin(), component.end(), '*');
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
