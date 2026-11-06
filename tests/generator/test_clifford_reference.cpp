// Garamon unit tests: self tests of the reference implementation tests/algebras/reference/CliffordReference.hpp
// (the generated libraries are checked against it, so it must be right)

#include <cstddef>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "CliffordReference.hpp"

using garamon_test::CliffordReference;
using Mv = CliffordReference::Mv;

namespace {

const std::vector<std::vector<double>> euclidean3 = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const std::vector<std::vector<double>> conformal3 = {
    {0, 0, 0, 0, -1}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {-1, 0, 0, 0, 0}};
const std::vector<std::vector<double>> nonOrthogonal3 = {{2, 1, 0}, {1, 1, 0}, {0, 0, -3}};
const std::vector<std::vector<double>> degenerate3 = {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}};

Mv randomMv(const CliffordReference &ref, std::mt19937 &rng) {
    std::uniform_int_distribution<int> coefficient(-2, 2);
    Mv mv = ref.zero();
    for(double &value : mv)
        value = coefficient(rng);
    return mv;
}

} // namespace

TEST_CASE("reference: Euclidean products of basis blades", "[reference]") {
    const CliffordReference ref(euclidean3);
    const unsigned int e1 = 1, e2 = 2, e3 = 4, e12 = 3, e13 = 5, e23 = 6, e123 = 7;

    CHECK(ref.bladeProduct(e1, e1) == ref.scalar(1.0));
    CHECK(ref.bladeProduct(e1, e2) == ref.blade(e12));
    CHECK(ref.bladeProduct(e2, e1) == ref.blade(e12, -1.0));
    CHECK(ref.bladeProduct(e12, e12) == ref.scalar(-1.0));
    CHECK(ref.bladeProduct(e123, e123) == ref.scalar(-1.0));
    CHECK(ref.bladeProduct(e12, e23) == ref.blade(e13));
    CHECK(ref.bladeProduct(e1, e23) == ref.blade(e123));
    CHECK(ref.bladeProduct(e3, e12) == ref.blade(e123));
    CHECK(ref.bladeProduct(e2, e13) == ref.blade(e123, -1.0));

    // e1 _| e123 = e23, e12 _| e123 = -e3
    CHECK(ref.leftContraction(ref.blade(e1), ref.blade(e123)) == ref.blade(e23));
    CHECK(ref.leftContraction(ref.blade(e12), ref.blade(e123)) == ref.blade(e3, -1.0));
    CHECK(ref.rightContraction(ref.blade(e123), ref.blade(e3)) == ref.blade(e12));
}

TEST_CASE("reference: non-orthogonal metric", "[reference]") {
    // e1.e2 = 1: e1 e2 = e1.e2 + e1^e2
    const CliffordReference ref(nonOrthogonal3);
    Mv expected = ref.blade(3);
    expected[0] = 1.0;
    CHECK(ref.bladeProduct(1, 2) == expected);
    expected[3] = -1.0;
    CHECK(ref.bladeProduct(2, 1) == expected);
    CHECK(ref.bladeProduct(1, 1) == ref.scalar(2.0));

    // (e1^e2)^2 = (e1.e2)^2 - e1^2 e2^2 = 1 - 2 = -1
    CHECK(ref.bladeProduct(3, 3) == ref.scalar(-1.0));
}

TEST_CASE("reference: conformal metric", "[reference]") {
    const CliffordReference ref(conformal3);
    const unsigned int e0 = 1, ei = 16;
    Mv expected = ref.blade(e0 | ei);
    expected[0] = -1.0;
    CHECK(ref.bladeProduct(e0, ei) == expected);
    CHECK(ref.bladeProduct(e0, e0) == ref.zero());
    CHECK(ref.bladeProduct(ei, ei) == ref.zero());
}

TEST_CASE("reference: algebraic identities", "[reference]") {
    for(const auto &metric : {euclidean3, conformal3, nonOrthogonal3, degenerate3}) {
        const CliffordReference ref(metric);
        std::mt19937 rng(1);
        for(unsigned int t = 0; t < 10; ++t) {
            const Mv a = randomMv(ref, rng), b = randomMv(ref, rng), c = randomMv(ref, rng);

            // associativity
            CHECK(CliffordReference::maxAbsDiff(ref.geometric(ref.geometric(a, b), c), ref.geometric(a, ref.geometric(b, c))) < 1.0e-9);
            CHECK(CliffordReference::maxAbsDiff(ref.outer(ref.outer(a, b), c), ref.outer(a, ref.outer(b, c))) < 1.0e-9);

            // reverse of a product
            CHECK(CliffordReference::maxAbsDiff(ref.reverse(ref.geometric(a, b)), ref.geometric(ref.reverse(b), ref.reverse(a))) < 1.0e-9);

            // vectors: v v = B(v,v)
            Mv v = ref.zero();
            for(unsigned int i = 0; i < ref.dimension(); ++i)
                v[std::size_t(1) << i] = a[std::size_t(1) << i];
            double expected = 0.0;
            for(unsigned int i = 0; i < ref.dimension(); ++i)
                for(unsigned int j = 0; j < ref.dimension(); ++j)
                    expected += v[std::size_t(1) << i] * v[std::size_t(1) << j] * metric[i][j];
            CHECK(CliffordReference::maxAbsDiff(ref.geometric(v, v), ref.scalar(expected)) < 1.0e-9);

            // the outer product does not depend on the metric
            const CliffordReference euclidean(std::vector<std::vector<double>>(metric.size(), std::vector<double>(metric.size(), 0.0)));
            CHECK(CliffordReference::maxAbsDiff(ref.outer(a, b), euclidean.outer(a, b)) < 1.0e-9);
        }
    }
}

TEST_CASE("reference: right complement", "[reference]") {
    const CliffordReference ref(degenerate3);
    for(unsigned int a = 0; a < 8; ++a)
        CHECK(ref.outer(ref.blade(a), ref.rightComplement(ref.blade(a))) == ref.pseudoScalar());
}
