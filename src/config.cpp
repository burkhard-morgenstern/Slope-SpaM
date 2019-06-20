#include "config.hpp"

#include <fmt/format.h>
#include <range/v3/action.hpp>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include "args.hpp"
#include "filesystem.hpp"

namespace fs = std::filesystem;
namespace ra = ranges::action;
namespace rv = ranges::view;

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

        if (results.size() == 0) {
            throw config_exception("No input files!\n");
        }

        return results;
    }

    auto parse_pattern(std::string const& patternflag)
        -> spam::pattern
    {
        auto const pattern = spam::pattern{patternflag};
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
        return wordlengths
            | rv::split(',')
            | rv::transform(
                [](auto&& rng) {
                    return std::string(&*rng.begin(),
                        ranges::distance(rng.begin(), rng.end()));
                })
            | rv::transform(
                [](auto&& length) {
                    return std::stoi(length);
                })
            | ranges::to<std::vector<size_t>>()
            | ra::sort
            | ra::unique;
    }

    config config::from_args(int argc, char** argv) {
        args::ArgumentParser parser("Slope-SpaM");
        args::ValueFlag<std::string> outfile(parser, "output file",
            "The output file. Ignored if multiple inputs are given.",
            {'o', "output"}, "");
        args::ValueFlag<std::string> patternflag(parser, "word pattern",
            "The binary word pattern used to create wordlists from sequences. "
            "May only include \'0\' and \'1\' characters.",
            {'p', "pattern"}, "111111111111111111111111111111111111");
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

        return {sanitize_inputs(input_files.Get()),
            outfile.Get(),
            parse_pattern(patternflag.Get()),
            as_reads.Get(),
            parse_wordlengths(wordlengths.Get())};
    }

}
