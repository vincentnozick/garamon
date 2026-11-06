// Garamon unit tests: m3ga with recursive products (tests/conf/rm3ga.conf)

#include <vector>

#include <rm3ga/Mvec.hpp> // IWYU pragma: keep

namespace ga = rm3ga;
#define GA_TAG "[rm3ga]"
#define GA_FULL_RANK 1
static const std::vector<std::vector<double>> gaMetric = {
    { 2, 0,  0},
    { 0, 1,  0},
    { 0, 0, -3}};

#include "AlgebraSuite.inc"
