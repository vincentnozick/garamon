// Garamon unit tests: src/MetricTools.cpp

#include <cmath>
#include <cstddef>
#include <vector>

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/Core>
#include <Eigen/LU> // IWYU pragma: keep (MatrixBase::determinant() is declared in Core but defined in LU)
#include <Eigen/SparseCore>

#include "MetricTools.hpp"
#include "Utilities.hpp"

namespace {

Eigen::MatrixXd c3gaMetric() {
    Eigen::MatrixXd m(5, 5);
    m << 0, 0, 0, 0, -1,
         0, 1, 0, 0, 0,
         0, 0, 1, 0, 0,
         0, 0, 0, 1, 0,
        -1, 0, 0, 0, 0;
    return m;
}

Eigen::MatrixXd p3gaMetric() {
    Eigen::MatrixXd m = Eigen::MatrixXd::Identity(4, 4);
    m(0, 0) = 0.0;
    return m;
}

Eigen::MatrixXd nonOrthogonalMetric() {
    Eigen::MatrixXd m(3, 3);
    m << 2, 1, 0,
         1, 1, 0,
         0, 0, -3;
    return m;
}

double maxAbs(const Eigen::MatrixXd &m) { return m.cwiseAbs().maxCoeff(); }

// determinant (minor) of the sub-matrix of 'm' made of the given rows and columns
double subdeterminant(const Eigen::MatrixXd &m, const std::vector<unsigned int> &rows, const std::vector<unsigned int> &cols) {
    Eigen::MatrixXd sub(rows.size(), cols.size());
    for(std::size_t i = 0; i < rows.size(); ++i)
        for(std::size_t j = 0; j < cols.size(); ++j)
            sub(i, j) = m(rows[i], cols[j]);
    return sub.determinant();
}

} // namespace

TEST_CASE("matrix predicates", "[metric_tools]") {
    const double epsilon = 1.0e-7;
    CHECK(isMatrixDiagonal(Eigen::MatrixXd::Identity(3, 3), epsilon));
    CHECK(isMatrixDiagonal(p3gaMetric(), epsilon));
    CHECK_FALSE(isMatrixDiagonal(c3gaMetric(), epsilon));

    CHECK(isMatrixIdentity(Eigen::MatrixXd::Identity(3, 3), epsilon));
    CHECK_FALSE(isMatrixIdentity(p3gaMetric(), epsilon));
    CHECK_FALSE(isMatrixIdentity(Eigen::MatrixXd::Identity(3, 4), epsilon));

    CHECK(isMatrixPermutationOfDiagonal(c3gaMetric(), epsilon));
    CHECK(isMatrixPermutationOfDiagonal(p3gaMetric(), epsilon));
    CHECK_FALSE(isMatrixPermutationOfDiagonal(nonOrthogonalMetric(), epsilon));
    Eigen::MatrixXd orthogonalColumns(2, 2); // M^T M is diagonal, but M is not a permutation of a diagonal matrix
    orthogonalColumns << 1, 1,
                         1, -1;
    CHECK_FALSE(isMatrixPermutationOfDiagonal(orthogonalColumns, epsilon));

    CHECK(getRank(c3gaMetric()) == 5u);
    CHECK(getRank(p3gaMetric()) == 3u);
    CHECK(getRank(nonOrthogonalMetric()) == 3u);
}

TEST_CASE("eigen decomposition of the metric", "[metric_tools]") {
    for(const Eigen::MatrixXd &metric : {c3gaMetric(), nonOrthogonalMetric(), p3gaMetric()}) {
        Eigen::MatrixXd P, D;
        eigenDecomposition(metric, P, D);
        CHECK(isMatrixDiagonal(D, 1.0e-12));
        CHECK(maxAbs(P * D * P.transpose() - metric) < 1.0e-9);

        // the refinement scales P and its inverse, the decomposition is unchanged
        Eigen::MatrixXd Pinv = P.transpose();
        const Eigen::MatrixXd scale = eigenRefinement(P, D, Pinv);
        CHECK(isMatrixDiagonal(scale, 1.0e-12));
        CHECK(maxAbs(P * D * Pinv - metric) < 1.0e-9);
        CHECK(maxAbs(Pinv * P - Eigen::MatrixXd::Identity(metric.rows(), metric.cols())) < 1.0e-9);
        CHECK(checkNumericalCleanUp(metric, P, D, Pinv, 1.0e-7));
        for(int i = 0; i < P.cols(); ++i)
            CHECK(minAbsNonZeroValue(P.col(i)) == Catch::Approx(1.0)); // smallest non-zero element of each column is 1
    }
}

TEST_CASE("minimum absolute non-zero value", "[metric_tools]") {
    Eigen::VectorXd v(4);
    v << 0.0, -0.5, 2.0, 0.0;
    CHECK(minAbsNonZeroValue(v) == 0.5);
    CHECK(minAbsNonZeroValue(Eigen::VectorXd::Zero(3)) == 0.0);
}

TEST_CASE("numerical clean up", "[metric_tools]") {
    const double epsilon = 1.0e-7;
    Eigen::VectorXd v(8);
    v << 1.0e-9,             // near zero
         0.9999999999,        // near integer
         -2.00000000001,      // near negative integer
         0.50000000001,       // near half
         2.0078125 + 1.0e-10, // near a multiple of 2^-7
         0.30000000001,       // near decimal
         0.123456789,         // no nice value: unchanged
         0.70710678118654757; // sqrt(0.5): unchanged

    const Eigen::VectorXd cleaned = vectorNumericalCleanUp(v, epsilon);
    CHECK(cleaned(0) == 0.0);
    CHECK(cleaned(1) == 1.0);
    CHECK(cleaned(2) == -2.0);
    CHECK(cleaned(3) == 0.5);
    CHECK(cleaned(4) == 2.0078125);
    CHECK(cleaned(5) == 0.3); // the nearest double
    CHECK(cleaned(6) == 0.123456789);
    CHECK(cleaned(7) == 0.70710678118654757);

    // matrix version: each row is cleaned
    Eigen::MatrixXd m(2, 8);
    m.row(0) = v.transpose();
    m.row(1) = -v.transpose();
    const Eigen::MatrixXd cleanedMatrix = numericalCleanUp(m, epsilon);
    CHECK(maxAbs(cleanedMatrix.row(0).transpose() - cleaned) == 0.0);
    CHECK(maxAbs(cleanedMatrix.row(1).transpose() + cleaned) == 0.0);
}

TEST_CASE("per grade transformation matrices are the compound matrices", "[metric_tools]") {
    const double epsilon = 1.0e-10;

    // identity
    for(unsigned int grade = 1; grade <= 4; ++grade) {
        const auto T = computePerGradeTransformationMatrix(Eigen::MatrixXd::Identity(4, 4), 4, grade, epsilon);
        CHECK(maxAbs(Eigen::MatrixXd(T) - Eigen::MatrixXd::Identity(bin_coeff(4, grade), bin_coeff(4, grade))) < 1.0e-12);
    }

    // any matrix: the element (l,m) of the grade k matrix is the minor (subdeterminant: rows of the l-th combination, columns of the m-th combination)
    Eigen::MatrixXd P(4, 4);
    P << 1, 2, 0, -1,
         0, 1, 3, 0,
         2, 0, 1, 1,
        -1, 1, 0, 2;
    for(unsigned int grade = 1; grade <= 4; ++grade) {
        INFO("grade " << grade);
        const auto T = computePerGradeTransformationMatrix(P, 4, grade, epsilon);
        const auto combinations = generateCombinations(4, grade);
        REQUIRE(T.rows() == (Eigen::Index)combinations.size());
        for(std::size_t l = 0; l < combinations.size(); ++l)
            for(std::size_t m = 0; m < combinations.size(); ++m)
                CHECK(T.coeff(l, m) == Catch::Approx(subdeterminant(P, combinations[l], combinations[m])).margin(1.0e-12));

        // no stored zero
        for(int k = 0; k < T.outerSize(); ++k)
            for(Eigen::SparseMatrix<double>::InnerIterator it(T, k); it; ++it)
                CHECK(std::fabs(it.value()) > epsilon);

        // inverse
        const auto Tinv = computeInverseTransformationMatrix(T, epsilon);
        CHECK(maxAbs(Eigen::MatrixXd(T * Tinv) - Eigen::MatrixXd::Identity(T.rows(), T.cols())) < 1.0e-9);
    }
}

TEST_CASE("all the transformation matrices of a metric", "[metric_tools]") {
    Eigen::MatrixXd P(3, 3);
    P << 1, 1, 0,
         1, -1, 0,
         0, 0, 1;

    std::vector<unsigned int> sizes, inverseSizes;
    std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> matrices, inverseMatrices;
    const auto components = computeTransformationMatricesToVector(P, 1.0e-7, sizes, inverseSizes, matrices, inverseMatrices);

    REQUIRE(matrices.size() == 4u);
    REQUIRE(inverseMatrices.size() == 4u);
    REQUIRE(sizes.size() == 4u);
    REQUIRE(inverseSizes.size() == 4u);

    unsigned int total = 0, inverseTotal = 0;
    for(unsigned int grade = 0; grade <= 3; ++grade) {
        INFO("grade " << grade);
        CHECK(matrices[grade].rows() == (Eigen::Index)bin_coeff(3, grade));
        CHECK(sizes[grade] == (unsigned int)matrices[grade].nonZeros());
        CHECK(inverseSizes[grade] == (unsigned int)inverseMatrices[grade].nonZeros());
        CHECK(maxAbs(Eigen::MatrixXd(matrices[grade] * inverseMatrices[grade]) - Eigen::MatrixXd::Identity(bin_coeff(3, grade), bin_coeff(3, grade))) < 1.0e-12);
        total += sizes[grade];
        inverseTotal += inverseSizes[grade];
    }

    // triplets (row, column, value)
    CHECK(components.first.size() == 3 * total);
    CHECK(components.second.size() == 3 * inverseTotal);
    CHECK(components.second[2 + 3 * 1] == Catch::Approx(0.5)); // first element of the grade 1 inverse
}
