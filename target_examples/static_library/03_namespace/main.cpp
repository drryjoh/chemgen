// main.cpp

#include "mylibrary.h"
#include <iostream>
using namespace std;
using namespace MyLibrary;

int main()
{
    // calling sayHello() function
    sayHello_();

    // calling addNumbers function and storing the result
    int result = MyLibrary::addNumbers_(5, 7);
    cout << "The result is: " << result << "\n";

    return 0;
}
