// Garamon unit tests: src/ConfigParser.cpp

#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include "ConfigParser.hpp"
#include "TestHelpers.hpp"

using garamon_test::TempDir;

TEST_CASE("configuration file loading", "[config_parser]") {
    TempDir dir("config_parser_file");

    SECTION("missing file") {
        CHECK_THROWS_AS(ConfigParser((dir.path() / "missing.conf").string()), std::runtime_error);
    }

    SECTION("Unix and Windows line endings") {
        const std::string content = "<namespace>\nc3ga\n</namespace>\n<metric>\n1 0\n0 -1\n</metric>\n<basis vector name>\n0 1\n</basis vector name>\n";
        std::string crlf;
        for(char c : content) {
            if(c == '\n')
                crlf += '\r';
            crlf += c;
        }
        for(const std::string &text : {content, crlf}) {
            const ConfigParser parser(dir.writeFile("test.conf", text, true));
            std::string name;
            REQUIRE(parser.readString("namespace", name));
            CHECK(name == "c3ga");
            Eigen::MatrixXd metric;
            REQUIRE(parser.readMatrix("metric", metric));
            CHECK(metric.rows() == 2);
            CHECK(metric.cols() == 2);
            CHECK(metric(1, 1) == -1.0);
            std::vector<std::string> names;
            REQUIRE(parser.readStringList("basis vector name", names));
            CHECK(names == std::vector<std::string>{"0", "1"});
        }
    }
}

TEST_CASE("tag extraction", "[config_parser]") {
    std::string value;

    SECTION("usual layout") {
        const auto parser = ConfigParser::fromString("# comment\n<name>\nvalue\n</name>\n");
        REQUIRE(parser.extract(value, "name"));
        CHECK(value == "value");
    }
    SECTION("same line, extra blanks") {
        const auto parser = ConfigParser::fromString("<name>  value \t</name>");
        REQUIRE(parser.extract(value, "name"));
        CHECK(value == "value");
    }
    SECTION("blank lines around the value") {
        const auto parser = ConfigParser::fromString("<name>\n\n  value\n\n</name>");
        REQUIRE(parser.extract(value, "name"));
        CHECK(value == "value");
    }
    SECTION("empty value") {
        const auto parser = ConfigParser::fromString("<name>\n</name>");
        REQUIRE(parser.extract(value, "name"));
        CHECK(value.empty());
    }
    SECTION("missing tag") {
        const auto parser = ConfigParser::fromString("<other>\nvalue\n</other>");
        CHECK_FALSE(parser.extract(value, "name"));
    }
    SECTION("missing closing tag") {
        const auto parser = ConfigParser::fromString("<name>\nvalue\n");
        CHECK_FALSE(parser.extract(value, "name"));
    }
    SECTION("tags whose name is a prefix of another tag") {
        const auto parser = ConfigParser::fromString("<metric decomposition refinement>\ntrue\n</metric decomposition refinement>\n<metric>\n1\n</metric>\n");
        REQUIRE(parser.extract(value, "metric"));
        CHECK(value == "1");
    }
}

TEST_CASE("unsigned integers", "[config_parser]") {
    unsigned int value = 42;
    CHECK(ConfigParser::fromString("<n>\n5\n</n>").readUInt("n", value));
    CHECK(value == 5u);
    CHECK(ConfigParser::fromString("<n> 17 </n>").readUInt("n", value));
    CHECK(value == 17u);

    CHECK_FALSE(ConfigParser::fromString("<n>\n-3\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("<n>\nabc\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("<n>\n5 6\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("<n>\n2.5\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("<n>\n\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("<n>\n99999999999999\n</n>").readUInt("n", value));
    CHECK_FALSE(ConfigParser::fromString("").readUInt("n", value));
}

TEST_CASE("floating point numbers", "[config_parser]") {
    double value = 0.0;
    CHECK(ConfigParser::fromString("<x>\n0.0000001\n</x>").readDouble("x", value));
    CHECK(value == 1.0e-7);
    CHECK(ConfigParser::fromString("<x>\n-2.5e3\n</x>").readDouble("x", value));
    CHECK(value == -2500.0);

    CHECK_FALSE(ConfigParser::fromString("<x>\nabc\n</x>").readDouble("x", value));
    CHECK_FALSE(ConfigParser::fromString("<x>\n0.5x\n</x>").readDouble("x", value));
    CHECK_FALSE(ConfigParser::fromString("<x>\n</x>").readDouble("x", value));
}

TEST_CASE("booleans", "[config_parser]") {
    bool value = false;
    CHECK(ConfigParser::fromString("<b>\ntrue\n</b>").readBool("b", value));
    CHECK(value);
    CHECK(ConfigParser::fromString("<b>\nFALSE\n</b>").readBool("b", value));
    CHECK_FALSE(value);
    CHECK(ConfigParser::fromString("<b>\nTrue\n</b>").readBool("b", value));
    CHECK(value);
    CHECK_FALSE(ConfigParser::fromString("<b>\nyes\n</b>").readBool("b", value));
    CHECK_FALSE(ConfigParser::fromString("<c>\ntrue\n</c>").readBool("b", value));
}

TEST_CASE("string lists", "[config_parser]") {
    std::vector<std::string> list;
    CHECK(ConfigParser::fromString("<l>\n0 1 2 3 i\n</l>").readStringList("l", list));
    CHECK(list == std::vector<std::string>{"0", "1", "2", "3", "i"});

    CHECK(ConfigParser::fromString("<l>\n  a   b\tc\n d  \n</l>").readStringList("l", list));
    CHECK(list == std::vector<std::string>{"a", "b", "c", "d"});

    CHECK(ConfigParser::fromString("<l>\n</l>").readStringList("l", list));
    CHECK(list.empty());

    CHECK_FALSE(ConfigParser::fromString("").readStringList("l", list));
}

TEST_CASE("matrices", "[config_parser]") {
    Eigen::MatrixXd m;

    SECTION("c3ga metric") {
        const auto parser = ConfigParser::fromString("<metric>\n 0  0  0  0 -1\n 0  1  0  0  0\n 0  0  1  0  0\n 0  0  0  1  0\n-1  0  0  0  0\n</metric>");
        REQUIRE(parser.readMatrix("metric", m));
        CHECK(m.rows() == 5);
        CHECK(m.cols() == 5);
        CHECK(m(0, 4) == -1.0);
        CHECK(m(4, 0) == -1.0);
        CHECK(m(2, 2) == 1.0);
    }
    SECTION("blank lines, tabs and decimals") {
        const auto parser = ConfigParser::fromString("<metric>\n\n1\t0.5\n\n0.5   -2.25  \n\n</metric>");
        REQUIRE(parser.readMatrix("metric", m));
        CHECK(m.rows() == 2);
        CHECK(m.cols() == 2);
        CHECK(m(0, 1) == 0.5);
        CHECK(m(1, 1) == -2.25);
    }
    SECTION("ragged matrix") {
        CHECK_FALSE(ConfigParser::fromString("<metric>\n1 0\n0\n</metric>").readMatrix("metric", m));
        CHECK_FALSE(ConfigParser::fromString("<metric>\n1\n0 1\n</metric>").readMatrix("metric", m));
    }
    SECTION("empty matrix") {
        CHECK_FALSE(ConfigParser::fromString("<metric>\n</metric>").readMatrix("metric", m));
        CHECK_FALSE(ConfigParser::fromString("<metric>\n\n\n</metric>").readMatrix("metric", m));
    }
    SECTION("not a number") {
        CHECK_FALSE(ConfigParser::fromString("<metric>\n1 x\n0 1\n</metric>").readMatrix("metric", m));
    }
    SECTION("missing matrix") {
        CHECK_FALSE(ConfigParser::fromString("").readMatrix("metric", m));
    }
}
