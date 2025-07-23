#include <iostream>
#include "MyUtility.h"

void MyUtility::sayHello() {
    std::cout << "Hello from MyUtility!" << std::endl;
}

int MyUtility::add(int a, int b) {
    return a + b;
}

double MyUtility::multiply(double x, double y) {
    return x * y;
}