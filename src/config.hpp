#include <filesystem>
#include <string>
#include <vector>

#include "spam/pattern.hpp"

namespace spam {

    class config_exception
        : public std::exception
    {
        std::string message;

    public:
        config_exception(std::string message)
            : message(std::move(message))
        {}

        auto what() const noexcept
            -> const char *
            override
        {
            return message.c_str();
        }
    };

    struct config {
        std::vector<std::filesystem::path> in;
        std::filesystem::path out;
        spam::pattern pattern;
        bool multi_fasta_as_reads;
        std::vector<size_t> wordlengths;

        static config from_args(int argc, char** argv);
    };

}
