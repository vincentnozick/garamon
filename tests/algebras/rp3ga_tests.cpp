// Garamon unit tests: p3ga with recursive products (tests/conf/rp3ga.conf)

#include <vector>

#include <rp3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = rp3ga;
#define GA_TAG "[rp3ga]"
#define GA_FULL_RANK 0
static const std::vector<std::vector<double>> gaMetric = {
    {0, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}};

#include "AlgebraSuite.inc"
