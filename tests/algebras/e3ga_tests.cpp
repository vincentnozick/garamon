// Garamon unit tests: Euclidean geometric algebra of R3 (conf/e3ga.conf)

#include <cmath>
#include <vector>

#include <e3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = e3ga;
#define GA_TAG "[e3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}};

#include "AlgebraSuite.inc"


TEST_CASE("e3ga: basis accessors", "[e3ga][accessors]") {
    CHECK(e3ga::E1 == 1u);
    CHECK(e3ga::E2 == 2u);
    CHECK(e3ga::E12 == 3u);
    CHECK(e3ga::E3 == 4u);
    CHECK(e3ga::E123 == 7u);

    const Mv e12 = e3ga::e12<double>();
    CHECK(e12.grades() == std::vector<unsigned int>{2});
    CHECK(e12[e3ga::E12] == 1.0);

    Mv mv;
    mv[e3ga::E1] = 2.0;
    mv[e3ga::E13] = 5.0;
    CHECK(mv.e1() == 2.0 * e3ga::e1<double>());
    CHECK(mv.e13() == 5.0 * e3ga::e13<double>());
    CHECK(mv.e2().isEmpty());
    CHECK(mv.e23().isEmpty());
}

TEST_CASE("e3ga: Euclidean geometry", "[e3ga][geometry]") {
    const Mv e1 = e3ga::e1<double>(), e2 = e3ga::e2<double>(), e3 = e3ga::e3<double>();
    const Mv I = e3ga::I<double>();

    SECTION("pseudo-scalar") {
        CHECK(sameMv(I * I, reference().scalar(-1.0)));
        CHECK(e3ga::Iinv<double>()[e3ga::E123] == -1.0);
    }

    SECTION("dual: A* = A _| I") {
        CHECK(sameMv(e1.dual(), toRef(e3ga::e23<double>())));
        CHECK(sameMv(e2.dual(), toRef(-e3ga::e13<double>())));
        CHECK(sameMv(e3ga::e12<double>().dual(), toRef(-e3)));
        CHECK(sameMv(Mv(1.0).dual(), toRef(I)));
    }

    SECTION("cross product a x b = -(a^b)*") {
        CHECK(sameMv(-(e1 ^ e2).dual(), toRef(e3)));
        CHECK(sameMv(-(e2 ^ e3).dual(), toRef(e1)));
        CHECK(sameMv(-(e3 ^ e1).dual(), toRef(e2)));
    }

    SECTION("rotor") {
        const double angle = 0.7;
        const Mv R = std::cos(angle / 2.0) - std::sin(angle / 2.0) * e3ga::e12<double>();
        CHECK(sameMv(R * ~R, reference().scalar(1.0)));
        const Mv rotated = R * e1 * ~R;
        Mv expected = std::cos(angle) * e1 + std::sin(angle) * e2;
        CHECK(maxAbsDiff(rotated, expected) < 1.0e-12);
        CHECK(maxAbsDiff(R * e3 * ~R, e3) < 1.0e-12);
    }

    SECTION("norm") {
        const Mv v = 3.0 * e1 + 4.0 * e2;
        CHECK(v.norm() == Catch::Approx(5.0));
        CHECK(v.quadraticNorm() == Catch::Approx(25.0));
    }
}
