// Garamon unit tests: conformal geometric algebra of R3 (conf/c3ga.conf)
// (see also c3ga_legacy_tests.cpp: hand-computed expectations, formerly plugin/c3gaUnitTests.cpp)

#include <vector>

#include <c3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = c3ga;
#define GA_TAG "[c3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 0, 0, 0, 0, -1},
    { 0, 1, 0, 0,  0},
    { 0, 0, 1, 0,  0},
    { 0, 0, 0, 1,  0},
    {-1, 0, 0, 0,  0}};

#include "AlgebraSuite.inc"

namespace {
    // conformal embedding of a Euclidean point
    Mv point(const double x, const double y, const double z) {
        return c3ga::e0<double>() + x * c3ga::e1<double>() + y * c3ga::e2<double>() + z * c3ga::e3<double>()
               + 0.5 * (x * x + y * y + z * z) * c3ga::ei<double>();
    }
}

TEST_CASE("c3ga: conformal geometry", "[c3ga][geometry]") {
    const Mv P = point(1, 2, 3), Q = point(-1, 0.5, 2);

    SECTION("points are null vectors and P.Q = -0.5 |p - q|^2") {
        CHECK(double(P | P) == Catch::Approx(0.0).margin(1.0e-12));
        CHECK(double(P | Q) == Catch::Approx(-0.5 * (4.0 + 2.25 + 1.0)));
    }

    SECTION("dual sphere and point membership") {
        // dual sphere of center c and radius r: c - 0.5 r^2 ei, a point is on the sphere iff P . S = 0
        const Mv center = point(1, 1, 1);
        const Mv dualSphere = center - 0.5 * 4.0 * c3ga::ei<double>(); // radius 2
        CHECK(double(point(3, 1, 1) | dualSphere) == Catch::Approx(0.0).margin(1.0e-12));
        CHECK(double(point(1, 1, -1) | dualSphere) == Catch::Approx(0.0).margin(1.0e-12));
        CHECK(double(point(1, 1, 1) | dualSphere) != Catch::Approx(0.0).margin(1.0e-12));
    }

    SECTION("a point lies on the sphere through 4 points") {
        const Mv sphere = point(1, 0, 0) ^ point(-1, 0, 0) ^ point(0, 1, 0) ^ point(0, 0, 1);
        CHECK(sphere.grades() == std::vector<unsigned int>{4});
        CHECK((point(0, -1, 0) ^ sphere).isEmpty());
        CHECK_FALSE((point(0, 0, 0) ^ sphere).isEmpty());
    }

    SECTION("translator") {
        // T = 1 - 0.5 t ei translates the points by t
        const Mv t = 1.0 * c3ga::e1<double>() + 2.0 * c3ga::e2<double>() - 1.0 * c3ga::e3<double>();
        const Mv T = 1.0 - 0.5 * t * c3ga::ei<double>();
        CHECK(maxAbsDiff(T * point(1, 2, 3) * ~T, point(2, 4, 2)) < 1.0e-12);
    }
}

TEST_CASE("c3ga: basis accessors", "[c3ga][accessors]") {
    CHECK(c3ga::E0 == 1u);
    CHECK(c3ga::E1 == 2u);
    CHECK(c3ga::Ei == 16u);
    CHECK(c3ga::E0123i == 31u);
    CHECK(c3ga::E12i == 22u);
}
