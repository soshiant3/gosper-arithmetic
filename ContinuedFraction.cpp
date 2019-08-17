//
// Created by Ina on 13.08.2019.
//

#include <cmath>
#include <iostream>
#include "ContinuedFraction.h"

//BivariateMoebiusTransform add {1,0,1,0,
//                               0,0,0,1};
//BivariateMoebiusTransform subtract {1,0,-1,0,
//                                    0,0,0,1};
//BivariateMoebiusTransform divide {1,0,0,0,
//                                  0,0,1,0};
//BivariateMoebiusTransform multiply {0,1,0,0,
//                                    0,0,0,1};

MoebiusTransform::MoebiusTransform(ContinuedFraction* x, int a, int b, int c, int d) : x(x), a(a), b(b), c(c), d(d) {}

const bool MoebiusTransform::must_feed() {
    if (c == 0 and d == 0) return false;
    if (c == 0 or d == 0) return true;
    else term = b / d;
    if (term == (int) (a / c)) {
        temp = b;
        b = d;
        d = temp - d * term;
        temp = a;
        a = c;
        c = temp - c * term;
        has_term = true;
        return false;
    }
    return true;
}
void MoebiusTransform::consume() {
    b = a;
    d = c;
}

void MoebiusTransform::consume(int next) {
    temp = b;
    b = a;
    a = temp + a * next;
    temp = d;
    d = c;
    c = temp + c * next;
}

const int MoebiusTransform::next() {
    has_term = false;
    return term;
}

const bool MoebiusTransform::has_next() {
    while (must_feed()) {
        if (x->has_next())
            consume(x->next());
        else consume();
    }
    return has_term;
}

// BIVARIATE

BivariateMoebiusTransform::BivariateMoebiusTransform(ContinuedFraction* x, ContinuedFraction* y, int a, int b, int c,
                                                     int d, int e, int f, int g, int h) : a(a), b(b), c(c), d(d), e(e),
                                                     f(f), g(g), h(h) {
    this->cfn[0] = x;
    this->cfn[1] = y;
}

int BivariateMoebiusTransform::choose_cfn() {
    return fabs(ae - dh) > fabs(cg - dh) ? 0 : 1;
}

const bool BivariateMoebiusTransform::must_feed() {
    if (e == 0 and f == 0 and g == 0 and h == 0) return false;
    if (h == 0) {
        cfn_index = g == 0 ? 0 : 1;
        return true;
    } else dh = ((double)d) / h;
    if (g == 0) {
        cfn_index = 1;
        return true;
    } else cg = ((double)c) / g;
    if (e == 0) {
        cfn_index = 0;
        return true;
    } else ae = ((double)a) / e;
    if (f == 0) {
        cfn_index = choose_cfn();
        return true;
    } else bf = ((double)b) / f;
    term = (int)dh;
    if (term == (int)ae and term == (int)cg and term == (int)bf) {
        temp = d;
        d = h;
        h = temp - h * term;
        temp = a;
        a = e;
        e = temp - e * term;
        temp = c;
        c = g;
        g = temp - g * term;
        temp = b;
        b = f;
        f = temp - f * term;
        has_term = true;
        return false;
    }
    cfn_index = choose_cfn();
    return true;
}

void BivariateMoebiusTransform::consume() {
    if (cfn_index == 0) {
        d = a;
        c = b;
        h = e;
        g = h;
    } else {
        d = c;
        a = b;
        h = g;
        e = f;
    }
}

void BivariateMoebiusTransform::consume(int next) {
    if (cfn_index == 0) {
        temp = d;
        d = a;
        a = temp + a * next;
        temp = c;
        c = b;
        b = temp + b * next;
        temp = h;
        h = e;
        e = temp + e * next;
        temp = g;
        g = f;
        f = temp + f * next;
    } else {
        temp = d;
        d = c;
        c = temp + c * next;
        temp = a;
        a = b;
        b = temp + b * next;
        temp = h;
        h = g;
        g = temp + g * next;
        temp = e;
        e = f;
        f = temp + f * next;
    }
}

const int BivariateMoebiusTransform::next() {
    has_term = false;
    return term;
}

const bool BivariateMoebiusTransform::has_next() {
    while (must_feed()) {
        if (cfn[cfn_index]->has_next())
            consume(cfn[cfn_index]->next());
        else consume();
    }
    return has_term;
}

const int ContinuedFraction::next() {}
const bool ContinuedFraction::has_next() {}
void Transform::consume() {}
const bool Transform::must_feed() {}
void Transform::consume(int) {}

ContinuedFraction* ContinuedFraction::operator+(ContinuedFraction &cf) {
    return new BivariateMoebiusTransform(this, &cf, 1, 0, 1, 0, 0, 0, 0, 1);
}
ContinuedFraction* ContinuedFraction::operator-(ContinuedFraction &cf) {
    return new BivariateMoebiusTransform(this, &cf, 1, 0, -1, 0, 0, 0, 0, 1);
}
ContinuedFraction* ContinuedFraction::operator*(ContinuedFraction &cf) {
    return new BivariateMoebiusTransform(this, &cf, 0, 1, 0, 0, 0, 0, 0, 1);
}
ContinuedFraction* ContinuedFraction::operator/(ContinuedFraction &cf) {
    return new BivariateMoebiusTransform(this, &cf, 1, 0, 0, 0, 0, 0, 1, 0);
}

std::ostream& operator<<(std::ostream &out, ContinuedFraction &cf) {
    int convergents[100];
    int i = 0;
    while (cf.has_next() && i < 100) {
        convergents[i++] = cf.next();
    }
    double val = convergents[--i];
    for (; i >= 0; --i) {
        val = convergents[i] + 1 / val;
    }
    return out << val;
}

Rational::Rational(int num, int denom) : num(num), denom(denom) {}

const bool Rational::has_next() {
    return denom != 0;
}

const int Rational::next() {
    int quotient = num / denom;
    int temp = num;
    num = denom;
    denom = temp % denom;
    return quotient;
}