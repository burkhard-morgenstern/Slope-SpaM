#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include <fmt/format.h>

#include "config.hpp"
#include "spam/distance_matrix.hpp"

namespace fs = std::filesystem;

class application {
  spam::config config;
  std::shared_ptr<ThreadPool> threadpool;

 public:
  explicit application(spam::config config)
      : config(std::move(config)),
        threadpool(std::make_shared<ThreadPool>(config.threads)) {}

 private:
  auto output_path(fs::path const& input_path) -> fs::path {
    return !fs::is_directory(input_path)
               ? fs::path{input_path}.replace_extension(".dmat")
               : fs::path{input_path.string() + ".dir.dmat"};
  }

  auto input_output_mapping() {
    auto result = std::vector<std::pair<fs::path, fs::path>>{};
    result.reserve(config.in.size());
    std::transform(config.in.begin(), config.in.end(),
                   std::back_inserter(result), [&](auto&& input_path) {
                     return std::make_pair(
                         input_path, (config.out == "" || config.in.size() != 1)
                                         ? output_path(input_path)
                                         : config.out);
                   });
    return result;
  }

 public:
  auto exec() -> int {
    if (!config.multi_fasta_as_reads) {
      for (auto const& [input_path, output_path] : input_output_mapping()) {
        fmt::print("Processing {} -> {}{}", input_path.string(),
                   output_path.string(), config.in.size() < 100 ? "\n" : "\r");
        auto sequences = load_sequences(input_path);
        if (sequences.size() == 0) {
          fmt::print(stderr, "Empty input path \"{}\".\n", input_path.string());
        } else {
          auto matrix = process(sequences);
          if (!matrix) {
            fmt::print("Couldn't calculate distance matrix!");
            return 1;
          }
          auto file = std::ofstream{output_path};
          file << *matrix;
          file.close();
        }
      }
    } else {
      fmt::print("Processing [");
      for (auto i = size_t{0}; i < config.in.size(); ++i) {
        fmt::print("{}", config.in[i].native());
        if (i + 1 < config.in.size()) {
          fmt::print(", ");
        }
      }
      fmt::print("] -> {}\n", config.out.string());
      auto sequences = std::vector<spam::sequence>{};
      std::transform(config.in.begin(), config.in.end(),
                     std::back_inserter(sequences), [&](auto&& input_path) {
                       return *spam::load_fasta_file(input_path);
                     });
      auto matrix = process(sequences);
      if (!matrix) {
        fmt::print("Couldn't calculate distance matrix!\n");
        return 1;
      }
      auto file = std::ofstream{config.out};
      file << *matrix;
      file.close();
    }

    return 0;
  }

 private:
  auto load_sequences(fs::path const& path) -> std::vector<spam::sequence> {
    if (fs::is_directory(path)) {
      return *spam::load_directory(path);
    } else {
      if (config.multi_fasta_as_reads) {
        return std::vector<spam::sequence>{*spam::load_fasta_file(path)};
      } else {
        auto assembled_sequences = *spam::load_multi_fasta_file(path);
        auto sequences = std::vector<spam::sequence>{};
        std::transform(assembled_sequences.begin(), assembled_sequences.end(),
                       std::back_inserter(sequences), [](auto&& seq) {
                         return spam::sequence{std::move(seq)};
                       });
        return sequences;
      }
    }
  }

  auto process(std::vector<spam::sequence> const& sequences)
      -> std::optional<spam::distance_matrix> {
    auto wordlengths = config.wordlengths;
    if (wordlengths.empty()) {
      if (spam::wordlist::kmax(sequences) <= spam::wordlist::kmin(sequences)) {
        return {};
      }
      wordlengths.push_back(spam::wordlist::kmin(sequences));
      wordlengths.push_back(spam::wordlist::kmax(sequences));
    }
    return spam::distance_matrix(std::move(sequences), config.pattern,
                                 wordlengths, threadpool);
  }
};

int main(int argc, char** argv) {

  
  // try {
  //   return application{spam::config::from_args(argc, argv)}.exec();
  // } catch (spam::config_exception const& e) {
  //   fmt::print(stderr, "{}", e.what());
  //   return 1;
  // } catch (spam::insufficient_pattern_exception const& e) {
  //   fmt::print(stderr, "Distance estimation failed:\n\t{}\n", e.what());
  // }
}
