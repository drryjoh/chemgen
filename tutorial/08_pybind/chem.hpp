#pragma once
#include <array>

constexpr int n_species = 5;
using Species = std::array<double, n_species>;

Species scale_species(const Species& y, double factor);

