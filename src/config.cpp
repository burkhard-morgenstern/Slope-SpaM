#include "config.hpp"

#include <fmt/format.h>
#include <range/v3/algorithm.hpp>

#include "args.hpp"
#include "filesystem.hpp"

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
            fmt::print("{}", parser.Help());
            exit(0);
        }
        catch (args::Error e)
        {
            fmt::print(stderr, "{}\n{}", e.what(), parser);
            exit(-1);
        }
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

        auto inputs = std::vector<fs::path>{};
        for (auto& in : input_files.Get()) {
            auto resolved = resolve_wildcards(in);
            std::copy(resolved.begin(), resolved.end(),
                std::back_inserter(inputs));
        }
        ranges::sort(inputs);

        return {inputs, outfile.Get(), patternflag.Get(), as_reads.Get()};
    }

}
