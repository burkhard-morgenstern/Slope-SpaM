#include "config.hpp"

#include <fstream>
#include <random>

#include <fmt/format.h>

#include "args.hpp"
#include "filesystem.hpp"
#include "spam/wordlist.hpp"
#include "string.hpp"

namespace fs = std::filesystem;

namespace spam {

    void parse_options(args::ArgumentParser & parser, int argc, char ** argv)
    {
        try
        {
            parser.ParseCLI(argc, argv);
        }
        catch (args::Help)
        {
            throw config_exception(fmt::format("{}", parser.Help()));
        }
        catch (args::Error e)
        {
            throw config_exception(fmt::format("{}\n", e.what()));
        }
    }

    auto sanitize_inputs(std::vector<std::string> inputs)
        -> std::vector<fs::path>
    {
        auto results = std::vector<fs::path>{};
        for (auto& in : inputs) {
            auto resolved = resolve_wildcards(in);
            std::copy(resolved.begin(), resolved.end(),
                std::back_inserter(results));
        }

        return results;
    }

    auto is_pattern(std::string const& str)
        -> bool
    {
        auto regex = std::regex(R"((0|1)*)");
        return regex_match(str, regex);
    }

    auto is_probability(std::string const& str)
        -> bool
    {
        auto regex = std::regex(R"(1\.0|0\.[0-9]+)");
        return regex_match(str, regex);
    }

    auto create_patternflag_with_coverage(double coverage)
        -> std::string
    {
        auto rd = std::random_device{};
        auto dist = std::uniform_real_distribution{0.0, 1.0};
        auto pattern = std::string{};
        for (auto match_count = size_t{0};
            match_count < spam::wordlist::max_wordsize();) {
            if (coverage > dist(rd)) {
                pattern += '1';
                ++match_count;
            } else {
                pattern += '0';
            }
        }
        return pattern;
    }

    auto random_pattern_from_flag(std::string const& patternflag)
        -> std::string
    {
        try {
            auto coverage = std::stod(patternflag);
            if (coverage < std::numeric_limits<double>::epsilon()) {
                throw config_exception(fmt::format(
                    "The given pattern coverage '{}' is not greater "
                    "than 0!", patternflag));
            }
            return create_patternflag_with_coverage(coverage);
        } catch (std::invalid_argument const&) {
            throw config_exception(fmt::format(
                "Conversion of the pattern '{}' to a converage value "
                "failed unexpectedly!", patternflag));
        } catch (std::out_of_range const&) {
            throw config_exception(fmt::format(
                "Pattern converage '{}' is not representable as a "
                "floating-point value!", patternflag));
        }
    }

    auto load_patternflag_from_file(fs::path const& filepath)
        -> std::string
    {
		auto file = std::ifstream{filepath};
        if (!file.is_open()) {
            throw config_exception(fmt::format(
                "Pattern file '{}' could not be opened!", filepath.native()));
        }

        std::string line;
        std::getline(file, line);
        if (!is_pattern(line)) {
            throw config_exception(fmt::format(
                "The first line of the pattern file '{}' is not a valid "
                "pattern!", filepath.native()));
        }
        return line;
    }

    auto resolve_patternflag_by_type(std::string const& patternflag)
        -> std::string
    {
        if (is_pattern(patternflag)) {
            return patternflag;
        }
        if (is_probability(patternflag)) {
            return random_pattern_from_flag(patternflag);
        }
        return load_patternflag_from_file(patternflag);
    }

    auto parse_pattern(std::string const& patternflag)
        -> spam::pattern
    {
        auto const pattern = spam::pattern{
            resolve_patternflag_by_type(patternflag)};
        if (pattern.weight() > spam::wordlist::max_wordsize()) {
            throw config_exception(fmt::format(
                "Unsupported pattern weight of {}! "
                "The supported maximum weight is {}!\n",
                pattern.weight(), spam::wordlist::max_wordsize()));
        }
        return pattern;
    }

    auto parse_wordlengths(std::string const& wordlengths)
        -> std::vector<size_t>
    {
        auto splitted = split(wordlengths, ',');
        auto result = std::vector<size_t>{};
        std::transform(splitted.begin(), splitted.end(),
            std::back_inserter(result),
            [](auto&& length) {
                return std::stoi(length);
            });
        std::sort(result.begin(), result.end());
        result.erase(
            std::unique(result.begin(), result.end()),
            result.end());
        return result;
    }

    auto config::from_args(int argc, char** argv)
        -> config
    {
        args::ArgumentParser parser("Slope-SpaM");
        args::ValueFlag<std::string> outfile(parser, "output file",
            "The output file. Ignored if multiple inputs are given.",
            {'o', "output"}, "");
        args::ValueFlag<std::string> patternflag(parser, "word pattern",
            "The binary word pattern used to create wordlists from sequences. "
            "May only include \'0\' and \'1\' characters.",
            {'p', "pattern"}, std::string(spam::wordlist::max_wordsize(), '1'));
        args::ValueFlag<std::string> wordlengths(parser, "word lengths",
            "Comma-separated list of wordlengths to consider when calculating"
            " the distance between two sequences.",
            {'k', "kmer-lengths"}, "");
        args::Flag as_reads(parser, "multi-fasta as reads",
            "Use sequences from multi-fasta files as reads of one sequence "
            "instead of multiple sequences.",
            {'a', "as-reads"}, false);
        args::HelpFlag help(parser, "help", "Show help.", {'h', "help"});
        args::PositionalList<std::string> input_files(parser, "input files",
            "Fasta files or directories of fasta files to process. If more than"
            " one file is given, the option output is ignored. Instead for each"
            " .fasta file a .dmat file and for each directory a .dir.dmat file"
            " with the same name is created.");
        parse_options(parser, argc, argv);

        auto inputs = sanitize_inputs(input_files.Get());
        std::sort(inputs.begin(), inputs.end());
        return {inputs,
            outfile.Get(),
            parse_pattern(patternflag.Get()),
            as_reads.Get(),
            parse_wordlengths(wordlengths.Get())};
    }

}
