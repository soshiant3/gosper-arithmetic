#include <iostream>
#include "ContinuedFraction.h"

int main() {
    Rational cf = Rational(141421, 100000);
    Rational z = Rational(5,2);
    ContinuedFraction* mt = cf / z;
    std::cout << *mt << std::endl;
    while(mt->has_next()) {
        std::cout << mt->next() << " ";
    }
    return 0;
}

template<int C, int R>
void gaussian(ContinuedFraction* mat[R][C]) {



}

