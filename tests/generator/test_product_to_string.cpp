// Garamon unit tests: src/ProductToString.cpp

#include <cstddef>
#include <list>
#include <string>
#include <vector>

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include "ConfigParser.hpp"
#include "MetaData.hpp"
#include "ProductToString.hpp"
#include "ProductTools.hpp"
#include "TestHelpers.hpp"
#include "Utilities.hpp"
#include "CliffordReference.hpp"

using namespace garamon_test;

namespace {

// all the numbers of a string, e.g. "{{ 1.5, -2}}" -> {1.5, -2}
std::vector<double> numbersIn(const std::string &text) {
    std::vector<double> numbers;
    std::string current;
    auto flush = [&]() {
        if(!current.empty() && current != "-" && current != ".")
            numbers.push_back(std::stod(current));
        current.clear();
    };
    for(char c : text) {
        if((c >= '0' && c <= '9') || c == '.' || c == '-' || c == 'e' || c == '+')
            current += c;
        else
            flush();
    }
    flush();
    return numbers;
}

MetaData metaDataFor(const std::string &metric, const std::string &names, const std::string &dimension) {
    return MetaData(ConfigParser::fromString(makeConfiguration("x", dimension, metric, names)));
}

} // namespace

TEST_CASE("constants to string", "[product_to_string]") {
    const MetaData c3ga(confFile("c3ga"));

    CHECK(basisVectorsToString(c3ga) == "{\"0\", \"1\", \"2\", \"3\", \"i\"}");
    CHECK(binomialCoefToString(5) == "binomialArray[6] = {1,5,10,10,5,1}");
    CHECK(reverseSignArrayToString(5) == "signReversePerGrade[6] = {1,1,-1,-1,1,1}");
    CHECK(reverseSignArrayToString(0) == "signReversePerGrade[1] = {1}");
    CHECK(perGradeStartingIndexToString({0, 1, 6, 16, 26, 31}) == "perGradeStartingIndex[6] = {0,1,6,16,26,31}");

    const std::string xorArrays = xorIndexToGradeAndHomogeneousIndexArraysToString(3, ProductTools(3));
    CHECK(xorArrays.find("xorIndexToGrade[] = {0,1,1,2,1,2,2,3}") != std::string::npos);
    CHECK(xorArrays.find("xorIndexToHomogeneousIndex[] = {0,0,1,0,2,1,2,0}") != std::string::npos);
}

TEST_CASE("metric to string", "[product_to_string]") {
    const MetaData c3ga(confFile("c3ga"));
    const std::string metric = metricToString(c3ga);
    CHECK(metric.find("e0\t0\t0\t0\t0\t-1\t") != std::string::npos);
    CHECK(metric.find("ei\t-1\t0\t0\t0\t0\t") != std::string::npos);

    // non-integer values are written without useless digits
    const MetaData m = metaDataFor("0.5 0\n0 -1.25", "a b", "2");
    const std::string decimal = metricToString(m);
    CHECK(decimal.find("ea\t0.5\t0\t") != std::string::npos);
    CHECK(decimal.find("eb\t0\t-1.25\t") != std::string::npos);
}

TEST_CASE("basis vector identifiers", "[product_to_string]") {
    const MetaData c3ga(confFile("c3ga"));
    CHECK(vectorIdToString(c3ga, false, 0) == "E0");
    CHECK(vectorIdToString(c3ga, false, 4) == "Ei");
    CHECK(vectorIdToString(c3ga, true, 0) == "e0");
    CHECK(vectorIdToString(c3ga, true, 4) == "ei");
}

TEST_CASE("component builder", "[product_to_string]") {
    const MetaData e3ga(confFile("e3ga"));
    const std::string constants = multivectorComponentBuilder(e3ga, constantsDefinition());
    CHECK(constants.find("const unsigned int E1 = 1;") != std::string::npos);
    CHECK(constants.find("const unsigned int E13 = 5;") != std::string::npos);
    CHECK(constants.find("const unsigned int E123 = 7;") != std::string::npos);
    CHECK(constants.find("project_") == std::string::npos);

    const std::string accessors = multivectorComponentBuilder(e3ga, oneComponentMultivectorPrototype());
    CHECK(accessors.find("Mvec e23() const {return this->extractOneComponent(2,3, 2);}") != std::string::npos);
}

TEST_CASE("right complement coefficients", "[product_to_string]") {
    // the generated library computes dual(e_A) = coefficients[k][p] e_B, with p = permutations[k][i] the position of e_B
    // (i: position of e_A in the k-vectors). For a degenerate metric, it must be the right complement: e_A ^ dual(e_A) = I
    for(unsigned int dimension = 3; dimension <= 6; ++dimension) {
        INFO("dimension " << dimension);
        const ProductTools tools(dimension);
        std::string components;
        const std::string code = fastRightComplementUtilities(dimension, tools, "", components);

        const std::size_t start = code.find("dualPermutations = ");
        REQUIRE(start != std::string::npos);
        const std::vector<double> permutations = numbersIn(code.substr(start + 19, code.find(";", start) - start - 19));
        const std::vector<double> coefficients = numbersIn(components);
        REQUIRE(permutations.size() == (std::size_t(1) << dimension));
        REQUIRE(coefficients.size() == (std::size_t(1) << dimension));

        const garamon_test::CliffordReference reference(std::vector<std::vector<double>>(dimension, std::vector<double>(dimension, 0.0)));
        unsigned int offset = 0;
        for(unsigned int grade = 0; grade <= dimension; ++grade) {
            for(unsigned int i = 0; i < bin_coeff(dimension, grade); ++i) {
                const unsigned int blade = tools.getXorIndex(grade, i);
                const unsigned int p = (unsigned int)permutations[offset + i];
                const unsigned int complement = tools.getXorIndex(dimension - grade, p);
                const double coefficient = coefficients[offset + p];
                INFO("blade " << blade << " -> " << coefficient << " * blade " << complement);
                CHECK(reference.outer(reference.blade(blade), reference.blade(complement, coefficient)) == reference.pseudoScalar());
            }
            offset += bin_coeff(dimension, grade);
        }
    }
}

TEST_CASE("coefficients are written with full precision", "[product_to_string]") {
    const double third = 1.0 / 3.0;

    SECTION("transformation matrices") {
        const std::vector<double> triplets = {0, 1, third, 1, 0, -third};
        const std::string text = perGradetransformMatricesToString(triplets, 1, 0, 2);
        const std::vector<double> numbers = numbersIn(text.substr(text.find('"')));
        REQUIRE(numbers.size() == 6u);
        CHECK(numbers[2] == third);
        CHECK(numbers[5] == -third);
    }

    SECTION("dual coefficients") {
        const std::string text = loadAllDualCoefficientsArray({0, 1}, doubleToString(third) + " " + doubleToString(-third) + " ");
        const std::size_t start = text.find("stringDualCoefficients=\"");
        REQUIRE(start != std::string::npos);
        const std::vector<double> numbers = numbersIn(text.substr(start + 24, text.find('"', start + 24) - start - 24));
        REQUIRE(numbers.size() == 2u);
        CHECK(numbers[0] == third);
    }

    SECTION("explicit products") {
        std::list<productComponent<double>> products = {{0, 0, 0, third}, {1, 1, 0, -2.0 * third}, {0, 1, 1, 1.0}};
        const std::string code = productListToString(products);
        CHECK(code.find(doubleToString(third) + "*mv1.coeff(0)*mv2.coeff(0)") != std::string::npos);
        CHECK(code.find(doubleToString(-2.0 * third) + "*mv1.coeff(1)*mv2.coeff(1)") != std::string::npos);
        CHECK(code.find("mv3.coeffRef(1) +=  mv1.coeff(0)*mv2.coeff(1);") != std::string::npos);
        CHECK(code.find("0.333333*") == std::string::npos);
    }

    SECTION("diagonal metric") {
        const MetaData m = metaDataFor("0.3 0\n0 -1", "a b", "2");
        const std::string text = diagonalMetricToString(m);
        CHECK(text.find("tmp<<0.29999999999999999,-1;") != std::string::npos);
    }
}

TEST_CASE("clean up of the product coefficients", "[product_to_string]") {
    // the basis changes (non-orthogonal metrics) leave small numerical errors in the coefficients
    std::list<productComponent<double>> products = {{0, 0, 0, 0.99999999999999967}, {0, 1, 0, -2.008412319093264e-16},
                                                    {1, 0, 0, -2.9999999999999996}, {1, 1, 1, 1.0 / 3.0}};
    cleanUpProductList(products, 1.0e-7);
    REQUIRE(products.size() == 3u);
    auto it = products.begin();
    CHECK(it->coefficient == 1.0);
    ++it;
    CHECK(it->coefficient == -3.0);
    ++it;
    CHECK(it->coefficient == 1.0 / 3.0); // no nice value close enough: unchanged
}
