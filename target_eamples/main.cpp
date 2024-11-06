#include <cmath>
#include "exp_gen.h"
#include "multiply_divide.h"
#include "pow_gen.h"
#include "constants.h"
#include "arrhenius.h"

#include <iostream>  // For printing the result to the console

int main() {
    // Call the arrhenius function with the specified parameters
    float result = arrhenius(float(100), float(1.3e6), float(1.5), float(1800));
    float dresult_dtemperature = darrhenius_dtemperature(float(100), float(1.3e6), float(1.5), float(1800));

    // Output the result
    std::cout << "Result of arrhenius(100, 1.3e6, 1.5, 1800): " << result << std::endl;
    std::cout << "Result of darrhenius_dtemperature(100, 1.3e6, 1.5, 1800): " << dresult_dtemperature << std::endl;

    return 0;
}
            