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
        args::ValueFlag<std::string> infile(parser, "Input file", "Specify input file name", {'i', "input"}, args::Options::Required);
        args::ValueFlag<std::string> outfile(parser, "Output file", "Specify output file name", {'o', "output"}, "out.dmat");
        args::ValueFlag<std::string> patternflag(parser, "Pattern", "Use this pattern if specified", {'p', "pattern"}, "111111111111111111111111111111111111");
        args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
        parse_options(parser, argc, argv);

        return {infile.Get(), outfile.Get(), patternflag.Get()};
    }

}