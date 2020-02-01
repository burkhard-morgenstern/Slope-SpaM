#include "distance_matrix.hpp"

#include <thread>

#include <fmt/format.h>

#include "../math.hpp"

namespace fs = std::filesystem;

namespace spam {

distance_matrix::distance_matrix(
    std::vector<spam::sequence> sequences, spam::pattern pattern,
    std::vector<size_t> wordlengths,
    std::optional<std::shared_ptr<ThreadPool>> threadpool /* = {} */)
    : sequences(std::move(sequences)),
      pattern(pattern),
      wordlengths(std::move(wordlengths)),
      threadpool(std::move(threadpool)) {
  calculate();
}

auto distance_matrix::size() const -> size_t { return sequences.size(); }

auto distance_matrix::column(size_t i) const
    -> std::pair<spam::sequence const&, std::vector<double> const&> {
  return {sequences[i], matrix[i]};
}

void distance_matrix::calculate() {
  initialize_matrix();
  create_wordlists();
  calculate_matrix();
}

void distance_matrix::initialize_matrix() {
  matrix.resize(sequences.size());
  for (auto& row : matrix) {
    row = std::vector<double>(sequences.size(), 0);
  }
}

void distance_matrix::create_wordlists() {
  auto max_wordlength =
      *std::max_element(wordlengths.begin(), wordlengths.end());
  try {
    pattern = pattern.reduce(max_wordlength);
  } catch (std::invalid_argument const& e) {
    throw insufficient_pattern_exception(
        fmt::format("The given pattern of weight {} is not sufficient for "
                    "the required maximal wordlength of {}!",
                    pattern.weight(), max_wordlength));
  }
  if (!threadpool) {
    create_wordlists_seq();
  } else {
    create_wordlists_par();
  }
}

void distance_matrix::create_wordlists_seq() {
  wordlists.reserve(sequences.size());
//   auto results = std::vector<std::future<spam::wordlist>>{};
  for (auto i = size_t{0}; i < sequences.size(); ++i) {
    wordlists.push_back(spam::wordlist(sequences[i], pattern));
  }
}

void distance_matrix::create_wordlists_par() {
  wordlists.reserve(sequences.size());
  auto results = std::vector<std::future<spam::wordlist>>{};
  for (auto i = size_t{0}; i < sequences.size(); ++i) {
    results.push_back((*threadpool)->enqueue([&, i = i]() {
      return spam::wordlist(sequences[i], pattern);
    }));
  }
  for (auto&& future : results) {
    wordlists.push_back(future.get());
  }
}

void distance_matrix::calculate_matrix() {
  if (!threadpool) {
    calculate_matrix_seq();
  } else {
    calculate_matrix_par();
  }
}

void distance_matrix::calculate_matrix_seq() {
  for (size_t i = 0; i < sequences.size() - 1; i++) {
    for (size_t j = i + 1; j < sequences.size(); j++) {
      auto const [probability, distance] = calculate_element(i, j);
      matrix[i][j] = distance;
      matrix[j][i] = distance;
    }
  }
}

void distance_matrix::calculate_matrix_par() {
  auto results = std::vector<std::future<std::pair<double, double>>>{};
  for (size_t i = 0; i < sequences.size() - 1; i++) {
    for (size_t j = i + 1; j < sequences.size(); j++) {
      results.push_back(
          (*threadpool)
              ->enqueue(
                  [&](size_t i, size_t j) { return calculate_element(i, j); },
                  i, j));
    }
  }

  size_t idx = 0;
  for (size_t i = 0; i < sequences.size() - 1; i++) {
    for (size_t j = i + 1; j < sequences.size(); j++) {
      auto const [probability, distance] = results[idx++].get();
      matrix[i][j] = distance;
      matrix[j][i] = distance;
    }
  }
}

auto jukes_cantor(double p) -> double {
  return -(3.0 / 4.0) * std::log(1.0 - (4.0 / 3.0) * (1.0 - p));
}

auto calculate_distance(std::vector<std::pair<size_t, size_t>> const& matches,
                        spam::sequence const& seq1, spam::sequence const& seq2)
    -> std::pair<double, double> {
  auto q = 0.25;
  auto values = std::vector<std::pair<long double, long double>>{};
  values.reserve(matches.size());
  for (auto& [k, count] : matches) {
    long double e =
        pow(q, k) * 2 * seq1.adjusted_size(k) * 2 * seq2.adjusted_size(k);
    values.emplace_back(k, log(count - e));
  }
  auto m = slope(values);
  auto p = exp(m) / ((1.0 - seq1.error_rate()) * (1.0 - seq2.error_rate()));
  // auto d = -(3.0 / 4.0) * log(1.0 - (4.0 / 3.0) * (1.0 - p));
  auto d = jukes_cantor(p);
  return {p, d};
}

auto distance_matrix::calculate_element(size_t i, size_t j) const
    -> std::pair<double, double> {
  auto matches = std::vector<std::pair<size_t, size_t>>{};
  matches.reserve(wordlengths.size());
  for (auto k : wordlengths) {
    matches.emplace_back(k, calculate_matches(wordlists[i], wordlists[j], k));
  }
  // auto kmax = *std::max_element(wordlengths.begin(), wordlengths.end());
  return calculate_distance(matches, sequences[i], sequences[j]);
}

std::ostream& operator<<(std::ostream& os, distance_matrix const& matrix) {
  os << matrix.size() << std::endl;
  for (int i = 0; i < matrix.size(); i++) {
    auto const& [sequence, distances] = matrix.column(i);
    auto length = sequence.name().length();
    if (length < 10) {
      os << sequence.name();
      for (auto i = length; length < 10; ++length) {
        os << ' ';
      }
    } else if (length > 10) {
      os << sequence.name().substr(0, 10);
    } else {
      os << sequence.name() << ' ';
    }
    for (int j = 0; j < matrix.size(); j++) {
      os << distances[j] << ' ';
    }
    os << "\n";
  }
  return os;
}

}  // namespace spam
