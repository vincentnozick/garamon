// Garamon unit tests: src/MetaData.cpp

#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <Eigen/Core>

#include "ConfigParser.hpp"
#include "MetaData.hpp"
#include "MetricTools.hpp"
#include "TestHelpers.hpp"

using namespace garamon_test;

namespace {

MetaData fromString(const std::string &content) { return MetaData(ConfigParser::fromString(content)); }

// the transformation matrix columns define an orthogonal basis: P^T M P = diag(diagonalMetric), P^-1 P = Id
void checkDiagonalization(const MetaData &metaData) {
    const unsigned int n = metaData.dimension;
    REQUIRE(metaData.diagonalMetric.size() == (Eigen::Index)n);
    if(metaData.inputMetricDiagonal) {
        CHECK((metaData.diagonalMetric - metaData.metric.diagonal()).cwiseAbs().maxCoeff() < 1.0e-12);
        return;
    }
    const Eigen::MatrixXd &P = metaData.transformationMatrix;
    const Eigen::MatrixXd &Pinv = metaData.inverseTransformationMatrix;
    REQUIRE(P.rows() == (Eigen::Index)n);
    const Eigen::MatrixXd orthogonalMetric = P.transpose() * metaData.metric * P;
    CHECK(isMatrixDiagonal(orthogonalMetric, 1.0e-6));
    CHECK((orthogonalMetric.diagonal() - metaData.diagonalMetric).cwiseAbs().maxCoeff() < 1.0e-6);
    CHECK((Pinv * P - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff() < 1.0e-6);
}

} // namespace

TEST_CASE("meta data of the tested algebras", "[metadata]") {

    SECTION("c3ga") {
        const MetaData metaData(confFile("c3ga"));
        CHECK(metaData.namespaceName == "c3ga");
        CHECK(metaData.dimension == 5u);
        CHECK(metaData.basisVectorName == std::vector<std::string>{"0", "1", "2", "3", "i"});
        CHECK_FALSE(metaData.inputMetricDiagonal);
        CHECK_FALSE(metaData.identityMetric);
        CHECK(metaData.fullRankMetric);
        CHECK(metaData.inputMetricPermutationOfDiagonal);
        CHECK(metaData.useEigenRefinement);
        CHECK(metaData.useNumericalCleanUp);
        CHECK(metaData.epsilon == 1.0e-7);
        CHECK(metaData.maxDimPrecomputedProducts == 256u);
        CHECK(metaData.maxDimBasisAccessor == 256u);
        checkDiagonalization(metaData);
    }

    SECTION("c5ga") {
        const MetaData metaData(confFile("c5ga"));
        CHECK(metaData.dimension == 7u);
        CHECK_FALSE(metaData.inputMetricDiagonal);
        CHECK(metaData.fullRankMetric);
        checkDiagonalization(metaData);
    }

    SECTION("e3ga") {
        const MetaData metaData(confFile("e3ga"));
        CHECK(metaData.dimension == 3u);
        CHECK(metaData.inputMetricDiagonal);
        CHECK(metaData.identityMetric);
        CHECK(metaData.fullRankMetric);
        CHECK_FALSE(metaData.inputMetricPermutationOfDiagonal); // initialized even though not computed for a diagonal metric
        checkDiagonalization(metaData);
    }

    SECTION("p3ga") {
        const MetaData metaData(confFile("p3ga"));
        CHECK(metaData.dimension == 4u);
        CHECK(metaData.inputMetricDiagonal);
        CHECK_FALSE(metaData.identityMetric);
        CHECK_FALSE(metaData.fullRankMetric);
        CHECK(metaData.diagonalMetric(0) == 0.0);
        checkDiagonalization(metaData);
    }

    SECTION("n3ga (non orthogonal)") {
        const MetaData metaData(testConfFile("n3ga"));
        CHECK_FALSE(metaData.inputMetricDiagonal);
        CHECK_FALSE(metaData.inputMetricPermutationOfDiagonal);
        CHECK(metaData.fullRankMetric);
        checkDiagonalization(metaData);
    }
}

TEST_CASE("all the configuration files of conf/ are valid", "[metadata]") {
    const std::string name = GENERATE(as<std::string>{},
        "c2ga", "c3ga", "c3ga2", "c3ga3", "c4ga", "c5ga", "c6ga", "c7ga", "c8ga", "c9ga", "c10ga", "cc2ga",
        "e2ga", "e3ga", "e4ga", "e5", "e7ga", "p2ga", "p3ga", "p3ga2", "qc2ga", "qc3ga", "st3ga", "st3ga2", "v3ga");
    INFO("configuration " << name);
    const MetaData metaData(confFile(name));
    CHECK(metaData.checkConsistency());
    CHECK(metaData.basisVectorName.size() == metaData.dimension);
    checkDiagonalization(metaData);
}

TEST_CASE("inconsistent configurations are rejected", "[metadata]") {
    const std::string identity3 = "1 0 0\n0 1 0\n0 0 1";

    // reference: a valid configuration
    CHECK_NOTHROW(fromString(makeConfiguration("e3ga", "3", identity3, "1 2 3")));

    SECTION("wrong number of basis vector names") {
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", identity3, "1 2")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", identity3, "1 2 3 4")), std::runtime_error);
    }
    SECTION("ambiguous basis vector names") {
        // e12 would be both "1","2" and "12"
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", identity3, "1 2 12")), std::runtime_error);
    }
    SECTION("basis vector names that are not part of a C++ identifier") {
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", identity3, "1 2 +")), std::runtime_error);
    }
    SECTION("metric dimension") {
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", "1 0\n0 1", "1 2 3")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", "1 0 0\n0 1 0", "1 2 3")), std::runtime_error);
    }
    SECTION("non symmetric metric") {
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "3", "1 2 0\n0 1 0\n0 0 1", "1 2 3")), std::runtime_error);
    }
    SECTION("dimension") {
        CHECK_THROWS_AS(fromString(makeConfiguration("e0ga", "0", "1", "")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "-3", identity3, "1 2 3")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("e3ga", "three", identity3, "1 2 3")), std::runtime_error);
    }
    SECTION("namespace") {
        CHECK_THROWS_AS(fromString(makeConfiguration("3ga", "3", identity3, "1 2 3")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("e3-ga", "3", identity3, "1 2 3")), std::runtime_error);
        CHECK_THROWS_AS(fromString(makeConfiguration("", "3", identity3, "1 2 3")), std::runtime_error);
    }
    SECTION("missing entries") {
        CHECK_THROWS_AS(fromString("<namespace>\ne3ga\n</namespace>\n"), std::runtime_error);
        CHECK_THROWS_AS(fromString(""), std::runtime_error);
    }
    SECTION("missing file") {
        CHECK_THROWS_AS(MetaData(confFile("this_algebra_does_not_exist")), std::runtime_error);
    }
}

TEST_CASE("metric diagonalization flags", "[metadata]") {
    SECTION("non-unit diagonal metric") {
        const MetaData metaData = fromString(makeConfiguration("m3ga", "3", "2 0 0\n0 1 0\n0 0 -3", "1 2 3"));
        CHECK(metaData.inputMetricDiagonal);
        CHECK_FALSE(metaData.identityMetric);
        CHECK(metaData.fullRankMetric);
    }
    SECTION("permutation of a diagonal metric") {
        const MetaData metaData = fromString(makeConfiguration("x2ga", "2", "0 1\n1 0", "1 2"));
        CHECK_FALSE(metaData.inputMetricDiagonal);
        CHECK(metaData.inputMetricPermutationOfDiagonal);
        CHECK(metaData.fullRankMetric);
        checkDiagonalization(metaData);
    }
    SECTION("degenerate non-diagonal metric") {
        const MetaData metaData = fromString(makeConfiguration("x3ga", "3", "1 1 0\n1 1 0\n0 0 1", "1 2 3"));
        CHECK_FALSE(metaData.fullRankMetric);
        checkDiagonalization(metaData);
    }
}
