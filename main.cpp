#include <iostream>
#include "math.h"
#include "ContinuedFraction.h"
#include <time.h>

template<int C, int R>
void gaussian (ContinuedFraction* (&mat)[R][C], ContinuedFraction* (&b)[R]) {

    int rank = 0;
    for (int j = 0; j < C; j++) {
        int pivot = -1;
        ContinuedFraction *curr = nullptr;
        for (int i = j; i < R; i++) {
            ContinuedFraction* abs = mat[i][j]->abs();
            if (!curr || (*abs != 0 && *curr > *abs)) {
                pivot = i;
                curr = abs;
            }
        }
        if (pivot >= 0) {
            ContinuedFraction* temp = b[pivot];
            b[pivot] = b[rank];
            b[rank] = temp;
            for (int i = j; i < C; i++) {
                temp = mat[pivot][i];
                mat[pivot][i] = mat[rank][i];
                mat[rank][i] = temp;
            }

            b[rank] = *b[rank] / *mat[rank][j];
            for (int i = C - 1; i >= j; i--) {
                mat[rank][i] = *mat[rank][i] / *mat[rank][j];
            }

            for (int k = 0; k < R; k++) {
                if (k == j) continue;
                b[k] = *b[k] - *(*mat[k][j] * *b[rank]);
                for (int i = C - 1; i >= j; i--) {
                    mat[k][i] = (*mat[k][i]) - *(*mat[k][j] * *mat[rank][i]);
                }
            }
            rank++;
        }
    }
}

template<int C, int R>
void gaussian (double (&mat)[R][C], double (&b)[R]) {

    int rank = 0;
    for (int j = 0; j < C; j++) {
        int pivot = -1;
        double curr = 0;
        for (int i = j; i < R; i++) {
            double abs_v = abs(mat[i][j]);
            if (!curr || (abs_v != 0 && curr > abs_v)) {
                pivot = i;
                curr = abs_v;
            }
        }
        if (pivot >= 0) {
            double temp = b[pivot];
            b[pivot] = b[rank];
            b[rank] = temp;
            for (int i = j; i < C; i++) {
                temp = mat[pivot][i];
                mat[pivot][i] = mat[rank][i];
                mat[rank][i] = temp;
            }

            b[rank] = b[rank] / mat[rank][j];
            for (int i = C - 1; i >= j; i--) {
                mat[rank][i] = mat[rank][i] / mat[rank][j];
            }

            for (int k = 0; k < R; k++) {
                if (k == j) continue;
                b[k] = b[k] - mat[k][j] * b[rank];
                for (int i = C - 1; i >= j; i--) {
                    mat[k][i] = mat[k][i] - mat[k][j] * mat[rank][i];
                }
            }
            rank++;
        }
    }
}

int main() {
    ContinuedFraction* e = new Factory([](int n) {
        if (n == 0) return 2;
        return (n+1) % 3 == 0 ? 2 * (n+1) / 3 : 1;
    });

    std::cout << *e << std::endl;

    ContinuedFraction* mat[4][4] = {
            {new Rational(1,1), new Rational(2,1), new Rational(-3,1), new Rational(-1,1)},
            {new Rational(0,1), new Rational(-3,1), new Rational(2,1), new Rational(6,1)},
            {new Rational(-3,1), new Rational(-1,1), new Rational(3,1), new Rational(1,1)},
            {new Rational(2,1), new Rational(3,1), new Rational(2,1), new Rational(-1,1)}
    };
    ContinuedFraction* b[4] = {
            new Rational(0,1),
            new Rational(-8,1),
            new Rational(0,1),
            new Rational(-8,1)
    };

    double mat_double[4][4] = {
            {1,2,-3,-1},
            {0,-3,2,6},
            {-3,-1,3,1},
            {2,3,2,-1}
    };
    double b_double[4] = {
            0,
            -8,
            0,
            -8
    };

    double start = clock();
    gaussian(mat, b);
    for (int i = 0; i < 4; i++) {
        std::cout << *b[i] << std::endl;
    }
    std::cout << "Program finished in " << (clock() - start) << "ms" << std::endl;
    start = clock();
    gaussian(mat_double, b_double);
    for (int i = 0; i < 4; i++) {
        std::cout << b_double[i] << std::endl;
    }
    std::cout << "Program finished in " << (clock() - start) << "ms" << std::endl;
    return 0;
}


