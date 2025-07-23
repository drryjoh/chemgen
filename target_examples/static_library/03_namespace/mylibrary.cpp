// mylibrary.cpp

#include "mylibrary.h"
#include <iostream>
using namespace std;

// function 1
void sayHello()
{
    cout << "Hello from the static library!\n";
}
// function 2
int addNumbers(int a, int b) { return a + b; }

void MyLibrary::sayHello_()
{
    sayHello();
}

int MyLibrary::addNumbers_(int a, int b)
{
    return
    100 +
    addNumbers(a, b);
}