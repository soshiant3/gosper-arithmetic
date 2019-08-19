#include <iostream>
#include "ContinuedFraction.h"

int main() {
    Rational cf = Rational(5, 2);
    Rational z = Rational(2, 1);
    ContinuedFraction* prod = cf * z;
    ContinuedFraction* sum = cf + z;
    ContinuedFraction* prod1 = prod->copy();
    ContinuedFraction* sum1 = sum->copy();
    while (sum1->has_next()) {
        std::cout << sum1->next() << " ";
    }
    std::cout << std::endl;
    while (prod1->has_next()) {
        std::cout << prod1->next() << " ";
    }
    std::cout << std::endl;

    std::cout << *sum << " " << *prod << std::endl;

    
    return 0;
}

template<int C, int R>
void gaussian (ContinuedFraction* mat[R][C], ContinuedFraction* b[R]) {

    for (int i = 0; i < C; i++) {
        int pivot = -1;
        for (int j = 0; j < R; j++) {
            if (*mat[j][i] != 0) {

            }
        }
    }

}

