// Garamon unit tests: src/Generator.cpp (generation of complete libraries)

#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "Generator.hpp"
#include "MetaData.hpp"
#include "TestHelpers.hpp"

using namespace garamon_test;

namespace {

// all the placeholders used in the template files (data/), e.g. project_namespace, cmake_project_name_original_case
const std::set<std::string> &templatePlaceholders() {
    static const std::set<std::string> placeholders = [] {
        std::set<std::string> result;
        const std::regex placeholder("(cmake_)?project_[A-Za-z0-9_]+");
        for(const auto &entry : fs::recursive_directory_iterator(sourceDir() / "data")) {
            if(!entry.is_regular_file())
                continue;
            const std::string content = readWholeFile(entry.path());
            for(auto it = std::sregex_iterator(content.begin(), content.end(), placeholder); it != std::sregex_iterator(); ++it)
                result.insert(it->str());
        }
        return result;
    }();
    return placeholders;
}

std::vector<std::string> expectedFiles(const std::string &name) {
    return {
        "CMakeLists.txt", "README.md", "LICENCE.txt", "cheatSheet.txt", "cmake/" + name + "Config.cmake.in",
        "doc/CMakeLists.txt", "doc/Doxyfile-html.cmake", "doc/HOWTO/HOWTO-build", "doc/HOWTO/HOWTO-doc", "doc/HOWTO/HOWTO-test",
        "sample/CMakeLists.txt", "sample/sample.py", "sample/src/main.cpp",
        "src/" + name + "/BasisTransformations.hpp", "src/" + name + "/Constants.hpp", "src/" + name + "/DualCoefficients.hpp",
        "src/" + name + "/Geometric.hpp", "src/" + name + "/GeometricExplicit.hpp", "src/" + name + "/Inner.hpp",
        "src/" + name + "/InnerExplicit.hpp", "src/" + name + "/Mvec.hpp", "src/" + name + "/Outer.hpp",
        "src/" + name + "/OuterExplicit.hpp", "src/" + name + "/PythonBindings.cpp", "src/" + name + "/Utility.hpp"};
}

} // namespace

TEST_CASE("the template placeholders are found", "[generation]") {
    const auto &placeholders = templatePlaceholders();
    CHECK(placeholders.count("project_namespace") == 1u);
    CHECK(placeholders.count("cmake_project_name_original_case") == 1u);
    CHECK(placeholders.size() > 20u);
}

TEST_CASE("generation of the tested algebras", "[generation]") {
    const std::string name = GENERATE(as<std::string>{}, "c3ga", "c5ga", "e3ga", "p3ga", "m3ga", "n3ga");
    const bool testAlgebra = (name == "m3ga" || name == "n3ga");
    const std::string conf = testAlgebra ? testConfFile(name) : confFile(name);
    INFO("algebra " << name);

    TempDir output("generation_" + name);
    const std::string library = generateLibrary(conf, templateDir(), output.str());
    CHECK(fs::equivalent(library, output.path() / ("garamon_" + name)));

    // all the files are generated
    for(const std::string &file : expectedFiles(name)) {
        INFO("file " << file);
        CHECK(fs::is_regular_file(fs::path(library) / file));
    }
    std::string upperCase = name;
    for(char &c : upperCase)
        c = (char)toupper(c);

    // header-only library, found by the sample with find_package(<name> CONFIG)
    CHECK_FALSE(fs::exists(fs::path(library) / "src" / name / "Mvec.cpp"));
    CHECK_FALSE(fs::exists(fs::path(library) / "sample" / "modules"));
    const std::string config = readWholeFile(fs::path(library) / "cmake" / (name + "Config.cmake.in"));
    CHECK(config.find("include(\"${CMAKE_CURRENT_LIST_DIR}/" + name + "Targets.cmake\")") != std::string::npos);
    CHECK(config.find("check_required_components(" + name + ")") != std::string::npos);

    // no placeholder is left in the generated files
    const auto &placeholders = templatePlaceholders();
    for(const auto &entry : fs::recursive_directory_iterator(library)) {
        if(!entry.is_regular_file())
            continue;
        const std::string content = readWholeFile(entry.path());
        for(const std::string &placeholder : placeholders) {
            if(content.find(placeholder) != std::string::npos) {
                INFO("file " << fs::relative(entry.path(), library).generic_string());
                FAIL_CHECK("placeholder '" << placeholder << "' not substituted");
            }
        }
    }

    // the library namespace is used
    const std::string mvec = readWholeFile(fs::path(library) / "src" / name / "Mvec.hpp");
    CHECK(mvec.find("namespace " + name + "{") != std::string::npos);
    CHECK(mvec.find("#ifndef " + upperCase + "_MULTI_VECTOR_HPP__") != std::string::npos);
}

TEST_CASE("generated constants", "[generation]") {
    TempDir output("generation_constants");

    SECTION("c3ga") {
        const std::string library = generateLibrary(confFile("c3ga"), templateDir(), output.str());
        const std::string constants = readWholeFile(fs::path(library) / "src" / "c3ga" / "Constants.hpp");
        CHECK(constants.find("pseudoScalarInverse = -1;") != std::string::npos);
        CHECK(constants.find("constexpr unsigned int algebraDimension = 5;") != std::string::npos);
        CHECK(constants.find("const unsigned int E0123i = 31;") != std::string::npos);
    }

    SECTION("e3ga: the dual is a left contraction onto I, the inverse pseudo-scalar is -I") {
        const std::string library = generateLibrary(confFile("e3ga"), templateDir(), output.str());
        const std::string constants = readWholeFile(fs::path(library) / "src" / "e3ga" / "Constants.hpp");
        CHECK(constants.find("pseudoScalarInverse = -1;") != std::string::npos);
    }

    SECTION("m3ga: metric diag(2,1,-3), I^-1 = I~ / det = -I / -6") {
        const std::string library = generateLibrary(testConfFile("m3ga"), templateDir(), output.str());
        const std::string constants = readWholeFile(fs::path(library) / "src" / "m3ga" / "Constants.hpp");
        CHECK(constants.find("pseudoScalarInverse = " + std::string("0.16666666666666666") + ";") != std::string::npos);
    }

    SECTION("fast dual availability") {
        const std::string c3ga = generateLibrary(confFile("c3ga"), templateDir(), output.str());
        CHECK(readWholeFile(fs::path(c3ga) / "src" / "c3ga" / "Constants.hpp").find("fastDualAvailable = true;") != std::string::npos);
        const std::string p3ga = generateLibrary(confFile("p3ga"), templateDir(), output.str());
        CHECK(readWholeFile(fs::path(p3ga) / "src" / "p3ga" / "Constants.hpp").find("fastDualAvailable = true;") != std::string::npos);
        // non-orthogonal metric that is not a permutation of a diagonal matrix: the dual of e1 is not a single blade
        const std::string n3ga = generateLibrary(testConfFile("n3ga"), templateDir(), output.str());
        CHECK(readWholeFile(fs::path(n3ga) / "src" / "n3ga" / "Constants.hpp").find("fastDualAvailable = false;") != std::string::npos);
    }

    SECTION("p3ga: degenerate metric, no inverse pseudo-scalar") {
        const std::string library = generateLibrary(confFile("p3ga"), templateDir(), output.str());
        const std::string constants = readWholeFile(fs::path(library) / "src" / "p3ga" / "Constants.hpp");
        CHECK(constants.find("pseudoScalarInverse") == std::string::npos);
        const std::string cheatSheet = readWholeFile(fs::path(library) / "cheatSheet.txt");
        CHECK(cheatSheet.find("right complement") != std::string::npos);
    }
}

TEST_CASE("generation errors", "[generation]") {
    TempDir output("generation_errors");

    SECTION("the library already exists") {
        generateLibrary(confFile("e3ga"), templateDir(), output.str());
        CHECK_THROWS_AS(generateLibrary(confFile("e3ga"), templateDir(), output.str()), std::runtime_error);
    }

    SECTION("invalid configuration: nothing is generated") {
        const std::string conf = output.writeFile("bad.conf", makeConfiguration("bad", "3", "1 0\n0 1", "1 2 3"));
        CHECK_THROWS_AS(generateLibrary(conf, templateDir(), output.str()), std::runtime_error);
        CHECK_FALSE(fs::exists(output.path() / "garamon_bad"));
    }

    SECTION("missing configuration file") {
        CHECK_THROWS_AS(generateLibrary(confFile("this_algebra_does_not_exist"), templateDir(), output.str()), std::runtime_error);
    }

    SECTION("missing template directory") {
        CHECK_THROWS_AS(generateLibrary(confFile("e3ga"), output.str() + "no_templates/", output.str()), std::runtime_error);
    }
}

TEST_CASE("generation without trailing separators in the directories", "[generation]") {
    TempDir output("generation_separators");
    std::string templates = templateDir();
    templates.pop_back(); // remove the trailing '/'
    std::string out = output.str();
    out.pop_back();

    const std::string library = generateLibrary(confFile("e3ga"), templates, out);
    CHECK(fs::is_regular_file(output.path() / "garamon_e3ga" / "src" / "e3ga" / "Mvec.hpp"));
    CHECK(fs::equivalent(library, output.path() / "garamon_e3ga"));
    CHECK_FALSE(fs::exists(fs::path(out + "garamon_e3ga")));
}

TEST_CASE("the generated explicit products are consistent", "[generation]") {
    // every function referenced by the arrays of function pointers is defined, and the coefficients are cleaned up
    // (c3ga2, c3ga3 and cc2ga are not generated here: their generation takes from 15 s to minutes)
    const std::string name = GENERATE(as<std::string>{},
        "c2ga", "c3ga", "c4ga", "c5ga", "e2ga", "e3ga", "e4ga", "e5", "p2ga", "p3ga", "p3ga2",
        "st3ga", "st3ga2", "v3ga");
    INFO("configuration " << name);
    TempDir output("generation_explicit_" + name);
    const std::string library = generateLibrary(confFile(name), templateDir(), output.str());
    const MetaData metaData(confFile(name));

    for(const std::string file : {"OuterExplicit.hpp", "InnerExplicit.hpp", "GeometricExplicit.hpp"}) {
        INFO("file " << file);
        const std::string content = readWholeFile(fs::path(library) / "src" / metaData.namespaceName / file);

        // references in the arrays of function pointers: "outer_1_2<T>", "geometric_1_2_1<T>", ...
        // (no std::regex: too slow on these large files)
        std::set<std::string> referenced;
        for(std::size_t pos = content.find("<T>"); pos != std::string::npos; pos = content.find("<T>", pos + 3)) {
            std::size_t start = pos;
            while(start > 0 && (std::isalnum((unsigned char)content[start - 1]) || content[start - 1] == '_'))
                --start;
            const std::string identifier = content.substr(start, pos - start);
            if(identifier.rfind("outer_", 0) == 0 || identifier.rfind("inner_", 0) == 0 || identifier.rfind("geometric_", 0) == 0)
                referenced.insert(identifier);
        }
        if(std::string(file) != "GeometricExplicit.hpp" || metaData.dimension > 2) // no "middle grade" geometric product in dimension 2
            CHECK_FALSE(referenced.empty());
        for(const std::string &function : referenced) {
            INFO("function " << function);
            CHECK(content.find("void " + function + "(") != std::string::npos);
        }

        // no numerical noise in the coefficients "c*mv1.coeff(i)*mv2.coeff(j)", such as 0.99999999999999967 or 2.0e-16
        unsigned int noisyCoefficients = 0;
        for(std::size_t pos = content.find("*mv1.coeff("); pos != std::string::npos; pos = content.find("*mv1.coeff(", pos + 1)) {
            std::size_t start = pos;
            while(start > 0 && (std::isdigit((unsigned char)content[start - 1]) || std::strchr(".e+-", content[start - 1]) != nullptr))
                --start;
            if(start == pos)
                continue;
            const double coefficient = std::stod(content.substr(start, pos - start));
            const double distanceToInteger = std::fabs(coefficient - std::round(coefficient));
            if(coefficient == 0.0 || (distanceToInteger > 0.0 && distanceToInteger < 1.0e-9))
                ++noisyCoefficients;
        }
        CHECK(noisyCoefficients == 0u);
    }
}
