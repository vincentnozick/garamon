// Garamon unit tests
// TestHelpers.hpp
//
// Licence MIT
//
// Common tools for the unit tests of the generator.

#ifndef GARAMON_TESTS_TEST_HELPERS_HPP
#define GARAMON_TESTS_TEST_HELPERS_HPP

#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <string>
#include <system_error>

namespace garamon_test {

namespace fs = std::filesystem;

/// root directory of the Garamon sources (conf/, data/, ...)
inline fs::path sourceDir() { return fs::path(GARAMON_SOURCE_DIR); }

/// configuration file of the conf/ directory, e.g. confFile("c3ga")
inline std::string confFile(const std::string &name) { return (sourceDir() / "conf" / (name + ".conf")).string(); }

/// configuration file of the tests/conf/ directory, e.g. testConfFile("m3ga")
inline std::string testConfFile(const std::string &name) { return (sourceDir() / "tests" / "conf" / (name + ".conf")).string(); }

/// template directory, with a trailing '/' as documented in the generator help
inline std::string templateDir() { return (sourceDir() / "data").generic_string() + "/"; }

/// a fresh empty directory, removed at destruction
class TempDir {
public:
    explicit TempDir(const std::string &name) : dir(fs::path(GARAMON_TEST_OUTPUT_DIR) / name) {
        fs::remove_all(dir);
        fs::create_directories(dir);
    }
    ~TempDir() {
        std::error_code ignored;
        fs::remove_all(dir, ignored);
    }
    TempDir(const TempDir &) = delete;
    TempDir &operator=(const TempDir &) = delete;

    const fs::path &path() const { return dir; }

    /// path of the directory as a string, with a trailing '/'
    std::string str() const { return dir.generic_string() + "/"; }

    /// write a text file in the directory, return its path
    std::string writeFile(const std::string &name, const std::string &content, const bool binary = false) const {
        const fs::path file = dir / name;
        std::ofstream stream(file, binary ? std::ios::binary : std::ios::out);
        stream << content;
        return file.string();
    }

private:
    fs::path dir;
};

inline std::string readWholeFile(const fs::path &file) {
    std::ifstream stream(file, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
}

/// content of a complete configuration file, with the given parts
inline std::string makeConfiguration(const std::string &namespaceName, const std::string &dimension, const std::string &metric,
                                     const std::string &basisVectorNames) {
    return "# test configuration\n"
           "<namespace>\n" + namespaceName + "\n</namespace>\n\n"
           "<dimension>\n" + dimension + "\n</dimension>\n\n"
           "<metric>\n" + metric + "\n</metric>\n\n"
           "<basis vector name>\n" + basisVectorNames + "\n</basis vector name>\n\n"
           "<metric decomposition refinement>\ntrue\n</metric decomposition refinement>\n\n"
           "<metric decomposition numerical cleanup>\ntrue\n</metric decomposition numerical cleanup>\n\n"
           "<metric decomposition numerical cleanup espilon>\n0.0000001\n</metric decomposition numerical cleanup espilon>\n\n"
           "<max dimension precomputed products>\n256\n</max dimension precomputed products>\n\n"
           "<max dimension basis accessor>\n256\n</max dimension basis accessor>\n";
}

} // namespace garamon_test

#endif // GARAMON_TESTS_TEST_HELPERS_HPP
