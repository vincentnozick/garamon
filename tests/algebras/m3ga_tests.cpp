// Garamon unit tests: diagonal non-normalized metric (tests/conf/m3ga.conf)

#include <vector>

#include <m3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = m3ga;
#define GA_TAG "[m3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 2, 0,  0},
    { 0, 1,  0},
    { 0, 0, -3}};

#include "AlgebraSuite.inc"


TEST_CASE("m3ga: inverse pseudo-scalar of a metric whose determinant is -6 or -3", "[m3ga][pseudoscalar]") {
    // I~ I = det(metric)  =>  I^{-1} = I~ / det(metric)
    const Mv I = m3ga::I<double>();
    const double determinant = double(~I * I);
    CHECK(determinant != Catch::Approx(1.0));
    CHECK(sameMv(m3ga::Iinv<double>(), toRef(~I * (1.0 / determinant))));
}
