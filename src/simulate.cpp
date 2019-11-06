#include "generator.hpp"
#include "spam/distance_matrix.hpp"

#include <fmt/format.h>

using namespace std::string_literals;

using csv_row = std::tuple<double, double, std::vector<double>>;

auto create_sequence_pair(size_t length, double distance)
    -> std::vector<spam::sequence>
{
    auto seq1 = DNA_SeqGenerator(length);
    auto seq2 = mutate(seq1, distance / 2.0);
    seq1 = mutate(seq1, distance / 2.0);
    return std::vector<spam::sequence>{
        spam::assembled_sequence{"", std::string(seq1.begin(), seq1.end())},
        spam::assembled_sequence{"", std::string(seq2.begin(), seq2.end())}
    };
}

void write_csv(std::vector<csv_row> rows, std::string filename) {
    auto file = std::ofstream(filename);
    for (auto& [key, expected, values] : rows) {
        file << key << " " << expected << " ";
        for (auto i = size_t{0}; i < values.size(); ++i) {
            file << values[i];
            if (i + 1 < values.size()) {
                file << " "s;
            }
        }
        file << "\n";
    }
}

auto const distances = std::vector<double>{
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1};
auto const sequence_length = size_t{10000};
auto const sample_size = size_t{20000};
auto const pattern = spam::pattern{"111111111111111111111111111111"};

void create_csv(std::string const& filename, std::vector<size_t> const& wordlengths) {
    auto results = std::vector<csv_row>{};
    for (auto const distance : distances) {
        fmt::print("Processing distance {}...\r", distance);
        std::flush(std::cout);
        auto row = std::vector<double>{};
        row.resize(sample_size);
        auto indices = std::vector<size_t>{};
        indices.reserve(sample_size);
        for (auto i = size_t{0}; i < sample_size; ++i) {
            indices.push_back(i);
        }
        std::transform(
            indices.begin(), indices.end(), row.begin(),
            [&](auto const i) {
                auto sequences = create_sequence_pair(sequence_length, distance);
                auto distance_matrix = spam::distance_matrix{
                    std::move(sequences), pattern, wordlengths};
                return distance_matrix.column(0).second[1];
            });
        results.push_back(std::make_tuple(distance, distance, std::move(row)));
    }
    write_csv(results, filename);
}

int main(int argc, char *argv[]) {
    create_csv("out.csv", {11, 14});
}
