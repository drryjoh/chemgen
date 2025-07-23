// mylibrary.h

#ifndef MYLIBRARY_H
#define MYLIBRARY_H
#include <iostream>
using namespace std;

void sayHello();
int addNumbers(int a, int b);

namespace MyLibrary
{
    void sayHello_();

    int addNumbers_(int a, int b);
}
#endif // MYLIBRARY_H