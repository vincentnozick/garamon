// Garamon unit tests: non-orthogonal metric (tests/conf/n3ga.conf)

#include <vector>

#include <n3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = n3ga;
#define GA_TAG "[n3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 2, 1,  0},
    { 1, 1,  0},
    { 0, 0, -3}};

#include "AlgebraSuite.inc"


TEST_CASE("n3ga: inverse pseudo-scalar of a metric whose determinant is -6 or -3", "[n3ga][pseudoscalar]") {
    // I~ I = det(metric)  =>  I^{-1} = I~ / det(metric)
    const Mv I = n3ga::I<double>();
    const double determinant = double(~I * I);
    CHECK(determinant != Catch::Approx(1.0));
    CHECK(sameMv(n3ga::Iinv<double>(), toRef(~I * (1.0 / determinant))));
}
