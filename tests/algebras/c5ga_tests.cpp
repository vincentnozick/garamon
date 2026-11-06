// Garamon unit tests: conformal geometric algebra of R5 (conf/c5ga.conf)

#include <array>
#include <vector>

#include <c5ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = c5ga;
#define GA_TAG "[c5ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 0, 0, 0, 0, 0, 0, -1},
    { 0, 1, 0, 0, 0, 0,  0},
    { 0, 0, 1, 0, 0, 0,  0},
    { 0, 0, 0, 1, 0, 0,  0},
    { 0, 0, 0, 0, 1, 0,  0},
    { 0, 0, 0, 0, 0, 1,  0},
    {-1, 0, 0, 0, 0, 0,  0}};

#include "AlgebraSuite.inc"

namespace {
    // conformal embedding of a Euclidean point
    Mv point(const std::array<double, 5> &x) {
        Mv p = c5ga::e0<double>();
        double squaredNorm = 0.0;
        const Mv e[5] = {c5ga::e1<double>(), c5ga::e2<double>(), c5ga::e3<double>(), c5ga::e4<double>(), c5ga::e5<double>()};
        for(unsigned int i = 0; i < 5; ++i) {
            p += x[i] * e[i];
            squaredNorm += x[i] * x[i];
        }
        return p + 0.5 * squaredNorm * c5ga::ei<double>();
    }
}

TEST_CASE("c5ga: conformal points", "[c5ga][geometry]") {
    const std::array<double, 5> x = {1, 2, -1, 0.5, 3};
    const std::array<double, 5> y = {0, -1, 2, 1, -2};
    const Mv P = point(x), Q = point(y);

    // points are null vectors
    CHECK(double(P | P) == Catch::Approx(0.0).margin(1.0e-12));
    CHECK(double(Q | Q) == Catch::Approx(0.0).margin(1.0e-12));

    // P . Q = -0.5 |x - y|^2
    double squaredDistance = 0.0;
    for(unsigned int i = 0; i < 5; ++i)
        squaredDistance += (x[i] - y[i]) * (x[i] - y[i]);
    CHECK(double(P | Q) == Catch::Approx(-0.5 * squaredDistance));

    // e0 . ei = -1, Minkowski plane
    CHECK(double(c5ga::e0<double>() | c5ga::ei<double>()) == -1.0);
    CHECK((c5ga::e0<double>() * c5ga::e0<double>()).isEmpty());
    CHECK((c5ga::ei<double>() * c5ga::ei<double>()).isEmpty());
}

TEST_CASE("c5ga: basis accessors", "[c5ga][accessors]") {
    CHECK(c5ga::E0 == 1u);
    CHECK(c5ga::Ei == 64u);
    CHECK(c5ga::E012345i == 127u);
    Mv mv;
    mv[c5ga::E05i] = 2.0;
    CHECK(mv.e05i() == 2.0 * c5ga::e05i<double>());
}
