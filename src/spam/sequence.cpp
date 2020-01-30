#include "sequence.hpp"

#include <algorithm>
#include <fstream>

namespace fs = std::filesystem;

namespace spam {

auto operator>>(std::istream& is, assembled_sequence& seq) -> std::istream& {
  std::getline(is, seq.name);
  std::string nucleotides;
  std::getline(is, nucleotides, '>');
  nucleotides.erase(std::remove_if(nucleotides.begin(), nucleotides.end(),
                                   [](const char c) { return isspace(c); }),
                    nucleotides.end());
  seq.nucleotides = std::move(nucleotides);
  return is;
}

auto sequence::name() const -> std::string {
  if (std::holds_alternative<assembled_sequence>(*this)) {
    return std::get<assembled_sequence>(*this).name;
  } else {
    return std::get<unassembled_sequence>(*this).name;
  }
}

auto sequence::size() const -> size_t {
  if (std::holds_alternative<assembled_sequence>(*this)) {
    return std::get<assembled_sequence>(*this).size();
  } else {
    auto& seq = std::get<unassembled_sequence>(*this);
    auto result = size_t{0};
    for (auto& read : seq.reads) {
      result += read.size();
    }
    return result;
  }
}

auto sequence::adjusted_size(size_t wordlength) const -> size_t {
  if (std::holds_alternative<assembled_sequence>(*this)) {
    return size() - wordlength + 1;
  } else {
    auto& seq = std::get<unassembled_sequence>(*this);
    return size() - seq.reads.size() * wordlength;
  }
}

auto sequence::error_rate() const -> long double {
  if (std::holds_alternative<assembled_sequence>(*this)) {
    return 0;
  } else {
    return std::get<unassembled_sequence>(*this).error_rate;
  }
}

auto count_nucleotide(spam::sequence const& sequence, char nucleotide)
    -> size_t {
  if (std::holds_alternative<assembled_sequence>(sequence)) {
    auto& nucleotides = std::get<assembled_sequence>(sequence).nucleotides;
    return std::count(nucleotides.begin(), nucleotides.end(), nucleotide);
  } else {
    auto& seq = std::get<unassembled_sequence>(sequence);
    auto result = size_t{0};
    for (auto& read : seq.reads) {
      result += std::count(read.begin(), read.end(), nucleotide);
    }
    return result;
  }
}

auto background_match_probability(spam::sequence const& seq1,
                                  spam::sequence const& seq2) -> double {
  auto result = 0.0;
  for (auto c : {'A', 'C', 'G', 'T'}) {
    result += 1.0 * count_nucleotide(seq1, c) / seq1.size() *
              (1.0 * count_nucleotide(seq2, c) / seq2.size());
  }
  return result;
}

auto load_directory(fs::path const& path)
    -> std::optional<std::vector<sequence>> {
  static std::vector<std::string> allowed_extensions = {".fasta", ".fa", ".fas",
                                                        ".fna", ".txt"};

  if (!fs::is_directory(path)) {
    return {};
  }

  auto result = std::vector<sequence>{};
  for (auto& entry : fs::directory_iterator(path)) {
    auto const has_allowed_extension =
        std::find(allowed_extensions.begin(), allowed_extensions.end(),
                  entry.path().extension()) != allowed_extensions.end();
    if (entry.is_regular_file() && has_allowed_extension) {
      result.push_back(*load_fasta_file(entry));
    }
  }
  return {result};
}

auto load_fasta_file(fs::path const& filename) -> std::optional<sequence> {
  auto file = std::ifstream(filename);
  if (!file.is_open()) {
    return {};
  }

  auto reads = std::vector<assembled_sequence>{};
  auto dummy = std::string{};
  std::getline(file, dummy, '>');
  while (!file.eof()) {
    auto seq = assembled_sequence{};
    file >> seq;
    reads.push_back(std::move(seq));
  }

  if (reads.size() == 1) {
    return {sequence{reads[0]}};
  } else {
    auto result = unassembled_sequence{};
    result.name = filename.stem();
    result.reads.resize(reads.size());
    result.error_rate = 0.0024; //todo ?
    for (auto i = size_t{0}; i < reads.size(); ++i) {
      result.reads[i] = std::move(reads[i].nucleotides);
    }
    return {sequence{result}};
  }
}

auto load_multi_fasta_file(fs::path const& filename)
    -> std::optional<std::vector<assembled_sequence>> {
  auto file = std::ifstream(filename);
  if (!file.is_open()) {
    return {};
  }
  
  auto result = std::vector<assembled_sequence>{};
  auto dummy = std::string{};
  std::getline(file, dummy, '>');
  while (!file.eof()) {
    auto seq = assembled_sequence{};
    file >> seq;
    result.push_back(std::move(seq));
  }
  return {result};
}

}  // namespace spam
