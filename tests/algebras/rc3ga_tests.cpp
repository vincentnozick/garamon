// Garamon unit tests: c3ga with recursive products (tests/conf/rc3ga.conf)

#include <vector>

#include <rc3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = rc3ga;
#define GA_TAG "[rc3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 0, 0, 0, 0, -1},
    { 0, 1, 0, 0,  0},
    { 0, 0, 1, 0,  0},
    { 0, 0, 0, 1,  0},
    {-1, 0, 0, 0,  0}};

#include "AlgebraSuite.inc"
