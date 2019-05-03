#include <string>

namespace spam {

    struct config {
        std::string in;
        std::string out;
        std::string pattern;

        static config from_args(int argc, char** argv);
    };

}
