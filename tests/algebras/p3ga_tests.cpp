// Garamon unit tests: projective geometric algebra of R3, degenerate metric (conf/p3ga.conf)

#include <vector>

#include <p3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = p3ga;
#define GA_TAG "[p3ga]"
#define GA_FULL_RANK 0
static const std::vector<std::vector<double>> gaMetric = {
    {0, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};

#include "AlgebraSuite.inc"


TEST_CASE("p3ga: degenerate metric", "[p3ga][metric]") {
    const Mv e0 = p3ga::e0<double>(), e1 = p3ga::e1<double>(), e2 = p3ga::e2<double>();

    CHECK((e0 * e0).isEmpty());
    CHECK((e0 | e0).isEmpty());
    CHECK((e0 | e1).isEmpty());
    CHECK(sameMv(e1 * e1, reference().scalar(1.0)));
    CHECK(sameMv(e0 * e1, toRef(p3ga::e01<double>())));
    CHECK(sameMv(e1 * e0, toRef(-p3ga::e01<double>())));

    // the pseudo-scalar is null, but the outer product is still defined
    const Mv I = p3ga::I<double>();
    CHECK((I * I).isEmpty());
    CHECK(sameMv(e0 ^ e1 ^ e2 ^ p3ga::e3<double>(), toRef(I)));
}

TEST_CASE("p3ga: right complement", "[p3ga][dual]") {
    // e_A ^ !e_A = I
    CHECK(sameMv(!p3ga::e0<double>(), toRef(p3ga::e123<double>())));
    CHECK(sameMv(!p3ga::e1<double>(), toRef(-p3ga::e023<double>())));
    CHECK(sameMv(!p3ga::e01<double>(), toRef(p3ga::e23<double>())));
    CHECK(sameMv(!Mv(1.0), toRef(p3ga::I<double>())));
    CHECK(sameMv(!p3ga::I<double>(), reference().scalar(1.0)));
}

TEST_CASE("p3ga: basis accessors", "[p3ga][accessors]") {
    CHECK(p3ga::E0 == 1u);
    CHECK(p3ga::E0123 == 15u);
    Mv mv;
    mv[p3ga::E02] = 3.0;
    CHECK(mv.e02() == 3.0 * p3ga::e02<double>());
    CHECK(mv.e01().isEmpty());
}
