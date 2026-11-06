// Garamon unit tests: src/Directory.cpp

#include <filesystem>
#include <stdexcept>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "Directory.hpp"
#include "TestHelpers.hpp"

using garamon_test::TempDir;

TEST_CASE("substitute", "[directory]") {
    std::string data = "namespace project_namespace{ project_namespace::Mvec }";
    substitute(data, "project_namespace", "c3ga");
    CHECK(data == "namespace c3ga{ c3ga::Mvec }");

    SECTION("pattern not found") {
        std::string text = "nothing to replace";
        substitute(text, "project_namespace", "c3ga");
        CHECK(text == "nothing to replace");
    }
    SECTION("the replacement is not interpreted") {
        std::string text = "a project_x b";
        substitute(text, "project_x", "$& $1 \\1 $$");
        CHECK(text == "a $& $1 \\1 $$ b");
    }
    SECTION("the pattern is not interpreted") {
        std::string text = "a.b axb";
        substitute(text, "a.b", "X");
        CHECK(text == "X axb");
    }
    SECTION("the replacement may contain the pattern") {
        std::string text = "ab ab";
        substitute(text, "ab", "abab");
        CHECK(text == "abab abab");
    }
    SECTION("non overlapping occurrences, from left to right") {
        std::string text = "aaa";
        substitute(text, "aa", "b");
        CHECK(text == "ba");
    }
    SECTION("empty pattern") {
        std::string text = "abc";
        substitute(text, "", "x");
        CHECK(text == "abc");
    }
}

TEST_CASE("path concatenation", "[directory]") {
    CHECK(joinPath("dir", "file") == "dir/file");
    CHECK(joinPath("dir/", "file") == "dir/file");
    CHECK(joinPath("dir\\", "file") == "dir\\file");
    CHECK(joinPath("", "file") == "file");
    CHECK(joinPath("a/b", "c/d") == "a/b/c/d");
}

TEST_CASE("files and directories", "[directory]") {
    TempDir dir("directory");

    SECTION("read and write") {
        const std::string file = dir.str() + "file.txt";
        const std::string content = "line 1\nline 2 \xc3\xa9\n";
        CHECK(writeFile(content, file));
        CHECK(readFile(file) == content);
        CHECK(directoryOrFileExists(file));
    }

    SECTION("read a missing file") {
        CHECK_THROWS_AS(readFile(dir.str() + "missing.txt"), std::runtime_error);
        CHECK_FALSE(directoryOrFileExists(dir.str() + "missing.txt"));
    }

    SECTION("write in a missing directory") {
        CHECK_FALSE(writeFile("x", dir.str() + "missing/file.txt"));
    }

    SECTION("make directories") {
        const std::string sub = dir.str() + "sub";
        CHECK_FALSE(directoryOrFileExists(sub));
        makeDirectory(sub);
        CHECK(garamon_test::fs::is_directory(sub));
        CHECK(directoryOrFileExists(sub));
        CHECK_THROWS_AS(makeDirectory(sub), std::runtime_error);               // already exists
        CHECK_THROWS_AS(makeDirectory(sub + "/a/b"), std::runtime_error);      // missing parent
    }
}
