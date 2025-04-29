#include "chem.hpp"

Species subtraction_from_species(const Species& y, double factor) {
    Species result = {};
    for (int i = 0; i < n_species; ++i) {
        result[i] = y[i] - factor;
    }
    return result;
}

Species scale_species(const Species& y, double factor) {
    Species result = {};
    for (int i = 0; i < n_species; ++i) {
        result[i] = y[i] * factor;
    }
    result  = subtraction_from_species(result, 1.0);
    return result;
}



