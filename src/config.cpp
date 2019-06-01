#include "config.hpp"

#include <iostream>

#include "args.hpp"

namespace spam {

    void parse_options(args::ArgumentParser & parser, int argc, char ** argv)
    {
        try
        {
            parser.ParseCLI(argc, argv);
        }
        catch (args::Help)
        {
            std::cout << parser;
            exit(0);
        }
        catch (args::Error e)
        {
            std::cerr << e.what() << std::endl;
            std::cerr << parser;
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
        args::HelpFlag help(parser, "help", "Show help.", {'h', "help"});
        args::PositionalList<std::string> input_files(parser, "input files",
            "Fasta files or directories of fasta files to process. If more than"
            " one file is given, the option output is ignored. Instead for each"
            " .fasta file a .dmat file and for each directory a .dir.dmat file"
            " with the same name is created.");
        parse_options(parser, argc, argv);

        return {infile.Get(), outfile.Get(), patternflag.Get()};
    }

}