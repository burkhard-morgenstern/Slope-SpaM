#include <filesystem>
#include <string>
#include <vector>

namespace spam {

    struct config {
        std::vector<std::filesystem::path> in;
        std::string out;
        std::string pattern;

        static config from_args(int argc, char** argv);
    };

}
