#include "wordlist.hpp"

#include <algorithm>
#include <cmath>

namespace fs = std::filesystem;

namespace spam {

auto wordlist::max_wordsize() -> size_t { return 4 * sizeof(word_t) - 1; }

auto wordlist::kmin(sequence const& seq) -> size_t {
  return std::ceil(std::log(seq.size()) / 0.87);
  // return std::ceil((std::log(seq.size()) + 0.69) / 0.875); //todo check
}

auto wordlist::kmin(std::vector<sequence> const& seqs) -> size_t {
  auto result = std::numeric_limits<size_t>::min();
  for (auto& seq : seqs) {
    result = std::max(result, kmin(seq));
  }
  return result;
}

auto wordlist::kmax(sequence const& seq) -> size_t {
  return std::floor(std::log(seq.size()) / 0.63);
  // return std::ceil(std::log(seq.size()) / 0.634); //todo check
}

auto wordlist::kmax(std::vector<sequence> const& seqs) -> size_t {
  auto result = std::numeric_limits<size_t>::max();
  for (auto& seq : seqs) {
    result = std::min(result, kmax(seq));
  }
  return result;
}

auto encode_sequence(std::string const& sequence) -> std::vector<int8_t> {
  auto result = std::vector<int8_t>{};
  result.reserve(sequence.size());
  std::transform(sequence.begin(), sequence.end(), std::back_inserter(result),
                 [](auto&& c) {
                   switch (toupper(c)) {
                     case 'A':
                       return 0;
                     case 'C':
                       return 1;
                     case 'T':
                       return 2;
                     case 'G':
                       return 3;
                    //  case 'C':
                    //    return 1;
                    //  case 'G':
                    //    return 2;
                    //  case 'T':
                    //    return 3;
                   }
                   return -1;
                 });
  return result;
}

template <class EncodedIt>
auto create_word(EncodedIt it, spam::pattern const& pattern) -> word_t {
  auto word = word_t{0};
  auto successful = true;
  for (auto pos : pattern) {
    if (it[pos] == -1) {
      successful = false;
      break;
    }
    word <<= 2;
    word += it[pos];
  }
  word <<= 8 * sizeof(word_t) - 2 * pattern.weight();
  return successful ? word : std::numeric_limits<word_t>::max();
//   return successful ? (word << 8 * sizeof(word_t) - 2 * pattern.weight())
//                     : std::numeric_limits<word_t>::max();
}

template <class EncodedIt>
auto create_revComp_word(EncodedIt it, spam::pattern const& pattern) -> word_t {
  auto revIt = it + pattern.size() - 1;

  auto word = word_t{0};
  auto successful = true;
  for (auto pos : pattern) {
    if (revIt[-pos] == -1) {
      successful = false;
      break;
    }
    word <<= 2;
    word += (revIt[-pos] + 2) % 4;
    // word += 3 - revIt[-pos];
  }
  word <<= 8 * sizeof(word_t) - 2 * pattern.weight();
  return successful ? word : std::numeric_limits<word_t>::max();
//   return successful ? (word << 8 * sizeof(word_t) - 2 * pattern.weight())
//                     : std::numeric_limits<word_t>::max();
}

auto create_words(std::string const& sequence, spam::pattern const& pattern)
    -> std::vector<word_t> {
  if (sequence.size() < pattern.size()) {
    return {};
  }
  
  auto encoded = encode_sequence(sequence);
  auto words = std::vector<word_t>{};
  words.reserve(2 * (sequence.size() - pattern.size() + 1));
  for (auto it = encoded.begin(); it != encoded.end() - pattern.size() + 1;
       ++it) {
    words.push_back(create_word(it, pattern));
    words.push_back(create_revComp_word(it, pattern));
  }
  words.erase(std::remove_if(words.begin(), words.end(),
                             [](auto const& word) {
                               return word ==
                                      std::numeric_limits<word_t>::max();
                             }),
              words.end());
  return words;
}

wordlist::wordlist(spam::sequence const& sequence, spam::pattern pattern)
    : pattern(pattern) {
  if (std::holds_alternative<assembled_sequence>(sequence)) {
    words = create_words(std::get<assembled_sequence>(sequence).nucleotides,
                         pattern);
  } else {
    auto& reads = std::get<unassembled_sequence>(sequence).reads;
    for (auto& read : reads) {
      auto tmp = create_words(read, pattern);
      std::copy(tmp.begin(), tmp.end(), std::back_inserter(words));
    }
  }
  std::sort(words.begin(), words.end());
}

auto calculate_matches(spam::wordlist const& wordlist1,
                       spam::wordlist const& wordlist2, size_t k) -> size_t {
  auto word_mask = std::numeric_limits<word_t>::max();
  word_mask <<= 8 * sizeof(word_t) - 2 * k;

  size_t count = 0;
  auto next_fn = [word_mask](auto it, auto&& wordlist) {
    return std::find_if_not(it, wordlist.end(), [it, word_mask](auto const& p) {
      return (p & word_mask) == ((*it) & word_mask);
    });
  };
  auto it1 = wordlist1.begin();
  auto next1 = next_fn(it1, wordlist1);
  auto it2 = wordlist2.begin();
  auto next2 = next_fn(it2, wordlist2);
  while (it1 != wordlist1.end() && it2 != wordlist2.end()) {
    auto const v1 = (*it1) & word_mask;
    auto const v2 = (*it2) & word_mask;
    if (v1 == v2) {
      count += std::distance(it1, next1) * std::distance(it2, next2);
    }
    if (v1 <= v2) {
      it1 = next1;
      next1 = next_fn(it1, wordlist1);
    }
    if (v1 >= v2) {
      it2 = next2;
      next2 = next_fn(it2, wordlist2);
    }
  }
  return count;
}

}  // namespace spam
