#include "MyUtility.h"
#include <iostream>

int main() {
    MyUtility::sayHello();
    std::cout << "3 + 4 = " << MyUtility::add(3, 4) << std::endl;
    std::cout << "2.5 * 4.0 = " << MyUtility::multiply(2.5, 4.0) << std::endl;
    return 0;
}