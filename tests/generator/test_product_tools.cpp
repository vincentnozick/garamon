// Garamon unit tests: src/ProductTools.cpp

#include <cmath>
#include <cstddef>
#include <functional>
#include <list>
#include <vector>

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Eigen/Core>

#include "ProductTools.hpp"
#include "Utilities.hpp"
#include "CliffordReference.hpp"

using garamon_test::CliffordReference;

namespace {

std::vector<std::vector<double>> diagonal(const std::vector<double> &values) {
    std::vector<std::vector<double>> metric(values.size(), std::vector<double>(values.size(), 0.0));
    for(std::size_t i = 0; i < values.size(); ++i)
        metric[i][i] = values[i];
    return metric;
}

Eigen::VectorXd toEigen(const std::vector<double> &values) {
    Eigen::VectorXd v(values.size());
    for(std::size_t i = 0; i < values.size(); ++i)
        v(i) = values[i];
    return v;
}

} // namespace

TEST_CASE("xor index <-> (grade, position) tables", "[product_tools]") {
    for(unsigned int dimension = 1; dimension <= 7; ++dimension) {
        const ProductTools tools(dimension);
        std::vector<unsigned int> countPerGrade(dimension + 1, 0);
        for(unsigned int xorIndex = 0; xorIndex < (1u << dimension); ++xorIndex) {
            INFO("dimension " << dimension << ", xor index " << xorIndex);
            const unsigned int grade = tools.getGrade(xorIndex);
            CHECK(grade == CliffordReference::grade(xorIndex));
            CHECK(tools.getHomogeneousIndex(xorIndex) < bin_coeff(dimension, grade));
            CHECK(tools.getXorIndex(grade, tools.getHomogeneousIndex(xorIndex)) == xorIndex);
            countPerGrade[grade]++;
        }
        for(unsigned int grade = 0; grade <= dimension; ++grade)
            CHECK(countPerGrade[grade] == bin_coeff(dimension, grade));
    }

    // order of the k-vector components (lexicographic order of the basis vector indices)
    const ProductTools tools(3);
    CHECK(tools.getXorIndex(2, 0) == 0b011u); // e12
    CHECK(tools.getXorIndex(2, 1) == 0b101u); // e13
    CHECK(tools.getXorIndex(2, 2) == 0b110u); // e23
    CHECK(tools.getXorIndex(3, 0) == 0b111u);
    CHECK(tools.getXorIndex(0, 0) == 0u);
}

TEST_CASE("blade product predicates and signs", "[product_tools]") {
    CHECK(ProductTools::hammingWeight(0xFFFFFFFFu) == 32u);

    for(unsigned int a = 0; a < 32; ++a)
        for(unsigned int b = 0; b < 32; ++b) {
            INFO("blades " << a << ", " << b);
            CHECK(ProductTools::hammingWeight(a) == CliffordReference::grade(a));
            CHECK(ProductTools::outerProductExists(a, b) == ((a & b) == 0));
            CHECK(ProductTools::productResultXorIndex(a, b) == (a ^ b));
            if((a & b) == 0)
                CHECK(ProductTools::outerProductSign(a, b) == (int)CliffordReference::outerSign(a, b));

            // the inner product of basis blades A.B is non null iff one is included in the other
            const unsigned int ga = CliffordReference::grade(a), gb = CliffordReference::grade(b);
            bool exists;
            if(ga < gb)      exists = (a & b) == a;
            else if(ga > gb) exists = (a & b) == b;
            else             exists = a == b;
            CHECK(ProductTools::innerProductExists(a, b) == exists);
        }
}

TEST_CASE("product coefficients with a diagonal metric", "[product_tools]") {
    const std::vector<double> metricValues = {2.0, 3.0, -5.0, 0.0};
    const CliffordReference reference(diagonal(metricValues));
    const Eigen::VectorXd metric = toEigen(metricValues);

    // with an orthogonal basis, e_A e_B = coefficient e_{A xor B}
    for(unsigned int a = 0; a < 16; ++a)
        for(unsigned int b = 0; b < 16; ++b) {
            INFO("e_" << a << " * e_" << b);
            const auto &product = reference.bladeProduct(a, b);
            CHECK(ProductTools::productCoefficientFromMetric(a, b, metric) == Catch::Approx(product[a ^ b]));
        }
}

TEST_CASE("explicit product lists with a diagonal metric", "[product_tools]") {
    const std::vector<double> metricValues = {2.0, -1.0, 3.0, 1.0};
    const unsigned int dimension = (unsigned int)metricValues.size();
    const CliffordReference reference(diagonal(metricValues));
    const Eigen::VectorXd metric = toEigen(metricValues);
    const ProductTools tools(dimension);

    for(unsigned int grade1 = 0; grade1 <= dimension; ++grade1)
        for(unsigned int grade2 = 0; grade2 <= dimension; ++grade2) {
            INFO("grades " << grade1 << ", " << grade2);

            // expected results: dense [index1][index2][xor index of the result]
            auto check = [&](const std::list<productComponent<double>> &list, const unsigned int gradeResult,
                             const std::function<CliffordReference::Mv(const CliffordReference::Mv &, const CliffordReference::Mv &)> &product) {
                std::vector<double> computed(bin_coeff(dimension, grade1) * bin_coeff(dimension, grade2) * (1u << dimension), 0.0);
                for(const auto &component : list) {
                    const unsigned int xorIndex3 = tools.getXorIndex(gradeResult, component.indexOfMv3);
                    computed[(component.indexOfMv1 * bin_coeff(dimension, grade2) + component.indexOfMv2) * (1u << dimension) + xorIndex3] += component.coefficient;
                }
                unsigned int failures = 0;
                for(unsigned int i1 = 0; i1 < bin_coeff(dimension, grade1); ++i1)
                    for(unsigned int i2 = 0; i2 < bin_coeff(dimension, grade2); ++i2) {
                        const auto expected = reference.gradePart(product(reference.blade(tools.getXorIndex(grade1, i1)),
                                                                          reference.blade(tools.getXorIndex(grade2, i2))), gradeResult);
                        for(unsigned int x = 0; x < (1u << dimension); ++x)
                            if(std::fabs(computed[(i1 * bin_coeff(dimension, grade2) + i2) * (1u << dimension) + x] - expected[x]) > 1.0e-12)
                                ++failures;
                    }
                CHECK(failures == 0u);
            };

            // outer product: grade1 + grade2
            if(grade1 + grade2 <= dimension) {
                const auto outer = tools.generateExplicitOuterProductList(grade1, grade2);
                CHECK(outer.size() == bin_coeff(dimension, grade1) * bin_coeff(dimension - grade1, grade2));
                check(outer, grade1 + grade2, [&reference](const CliffordReference::Mv &x, const CliffordReference::Mv &y) { return reference.outer(x, y); });
            }

            // inner product: |grade1 - grade2|
            const unsigned int gradeInner = grade1 > grade2 ? grade1 - grade2 : grade2 - grade1;
            const auto inner = tools.generateExplicitInnerProductListEuclideanSpace(grade1, grade2, metric);
            check(inner, gradeInner, [&reference](const CliffordReference::Mv &x, const CliffordReference::Mv &y) { return reference.dot(x, y); });

            // geometric product: the other grades
            const auto geometric = tools.generateExplicitGeometricProductListEuclideanSpace(grade1, grade2, metric);
            REQUIRE(geometric.size() == dimension + 1);
            for(unsigned int gradeResult = 0; gradeResult <= dimension; ++gradeResult) {
                if(gradeResult == grade1 + grade2 || gradeResult == gradeInner) {
                    CHECK(geometric[gradeResult].empty());
                    continue;
                }
                check(geometric[gradeResult], gradeResult, [&reference](const CliffordReference::Mv &x, const CliffordReference::Mv &y) { return reference.geometric(x, y); });
            }
        }
}

TEST_CASE("scale of the inverse pseudo-scalar", "[product_tools]") {
    Eigen::MatrixXd c3ga(5, 5);
    c3ga << 0, 0, 0, 0, -1,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
           -1, 0, 0, 0, 0;
    CHECK(getScaleInversePseudoScalar(c3ga) == Catch::Approx(-1.0));

    CHECK(getScaleInversePseudoScalar(Eigen::MatrixXd::Identity(3, 3)) == Catch::Approx(1.0));

    // the scale is 1/det(metric), not det(metric)
    Eigen::MatrixXd scaled = Eigen::MatrixXd::Zero(3, 3);
    scaled.diagonal() << 2, 1, -3;
    CHECK(getScaleInversePseudoScalar(scaled) == Catch::Approx(-1.0 / 6.0));

    Eigen::MatrixXd nonOrthogonal(3, 3);
    nonOrthogonal << 2, 1, 0,
                     1, 1, 0,
                     0, 0, -3;
    CHECK(getScaleInversePseudoScalar(nonOrthogonal) == Catch::Approx(-1.0 / 3.0));

    // degenerate metric
    Eigen::MatrixXd p3ga = Eigen::MatrixXd::Identity(4, 4);
    p3ga(0, 0) = 0;
    CHECK(getScaleInversePseudoScalar(p3ga) == 0.0);
}
