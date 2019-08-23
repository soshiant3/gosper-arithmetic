//
// Created by Ina on 13.08.2019.
//

#include <cmath>
#include <iostream>
#include "ContinuedFraction.h"
#include "IOConfig.cpp"

//BivariateMoebiusTransform add {1,0,1,0,
//                               0,0,0,1};
//BivariateMoebiusTransform subtract {1,0,-1,0,
//                                    0,0,0,1};
//BivariateMoebiusTransform divide {1,0,0,0,
//                                  0,0,1,0};
//BivariateMoebiusTransform multiply {0,1,0,0,
//                                    0,0,0,1};

MoebiusTransform::MoebiusTransform(ContinuedFraction* x, int a, int b, int c, int d) : x(x->copy()), a(a), b(b), c(c), d(d) {}

MoebiusTransform::~MoebiusTransform() {
    delete x;
}

const bool MoebiusTransform::must_feed() {
    if (c == 0 and d == 0) return false;
    if (c == 0 or d == 0) return true;
    else term = b / d;
    if (term == a / c) {
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

MoebiusTransform* MoebiusTransform::copy() const {
    return new MoebiusTransform(x->copy(), a, b, c, d);
}

// BIVARIATE

BivariateMoebiusTransform::BivariateMoebiusTransform(ContinuedFraction* x, ContinuedFraction* y, int a, int b, int c,
                                                     int d, int e, int f, int g, int h) : a(a), b(b), c(c), d(d), e(e),
                                                     f(f), g(g), h(h) {
    this->cfn[0] = x->copy();
    this->cfn[1] = y->copy();
}

BivariateMoebiusTransform::~BivariateMoebiusTransform() {
    delete cfn[0];
    delete cfn[1];
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
        g = f;
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

BivariateMoebiusTransform* BivariateMoebiusTransform::copy() const {
    return new BivariateMoebiusTransform(cfn[0], cfn[1], a, b, c, d, e, f, g, h);
}

const int ContinuedFraction::next() {}
const bool ContinuedFraction::has_next() {}
ContinuedFraction* ContinuedFraction::copy() const {}
void Transform::consume() {}
const bool Transform::must_feed() {}
void Transform::consume(int) {}

ContinuedFraction* ContinuedFraction::abs() {
    ContinuedFraction* a = copy();
    if (!a->has_next()) return a;
    return *this < 0 ? (delete a, new MoebiusTransform(this, -1, 0, 0, 1)) : a;
}

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
bool ContinuedFraction::operator==(int i) {
    ContinuedFraction* c = copy();
    if (c->has_next()) {
        int next = c->next();
        bool has_next = c->has_next();
        delete c;
        return i == next && !has_next;
    }
    delete c;
    return false; // this cf resembles infinity.
}
bool ContinuedFraction::operator!=(int i) {
    return ! (*this == i);
}

bool ContinuedFraction::operator==(ContinuedFraction &cf) {
    ContinuedFraction* a = copy();
    ContinuedFraction* b = cf.copy();
    for (int i = 0; i < iterations; i++) {
        if (a->has_next() && b->has_next()) {
            if (a->next() != b->next()) {
                delete a;
                delete b;
                return false;
            }
        } else {
            bool a_has_next = a->has_next(), b_has_next = b->has_next();
            delete a;
            delete b;
            return !a_has_next && !b_has_next;
        }
    }
    delete a;
    delete b;
    return true;
}

bool ContinuedFraction::operator!=(ContinuedFraction &cf) {
    return ! (*this == cf);
}

bool ContinuedFraction::operator<(ContinuedFraction &cf) {
    ContinuedFraction* a = copy();
    ContinuedFraction* b = cf.copy();
    for (int i = 0; i < iterations; i++) {
        if (a->has_next() && b->has_next()) {
            int a_next = a->next();
            int b_next = b->next();
            if (a_next < b_next) {
                delete a;
                delete b;
                return i % 2 == 0;
            } else if (a_next > b_next) {
                delete a;
                delete b;
                return i % 2 != 0;
            }
        } else {
            bool has_next = a->has_next();
            delete a;
            delete b;
            return has_next ? i % 2 == 0 : i % 2 != 0;
        }
    }
    delete a;
    delete b;
    return false;
}

bool ContinuedFraction::operator>(ContinuedFraction &cf) {
    return cf < *this;
}

bool ContinuedFraction::operator<(int i) {
    ContinuedFraction* a = copy();
    if (!a->has_next()) {
        delete a;
        return false;
    }
    int next = a->next();
    delete a;
    return next < i;
}

bool ContinuedFraction::operator>(int i) {
    ContinuedFraction* a = copy();
    if (!a->has_next()) {
        delete a;
        return true;
    }
    int a_next = a->next();
    bool has_next = a->has_next();
    delete a;
    return a_next > i || (has_next && a_next == i);
}

std::ostream& operator<<(std::ostream &out, ContinuedFraction &cf) {
    ContinuedFraction* cf0 = cf.copy();
    int P = 1, p = 0;
    int Q = 0, q = 1;
    int n = -1;
    unsigned int err = 0;
    bool komma = true;
    bool step = false;
    while (cf0->has_next() && (err < machine_eps || !step)) {
        if (err >= machine_eps) step = true;
        int next = cf0->next(), p_n = p, q_n = q;
        p = P;
        P = next * P + p_n;
        q = Q;
        Q = next * Q + q_n;
        err = Q; // sqrt of the error.
        n++;
    }
    if (!cf0->has_next()) {
        delete cf0;
        return Q == 1 ? out << P : (Q == 0 ? out << "inf" : out << P << "/" << Q);
    }
    if (!decimal) {
        delete cf0;
        return n % 2 == 0 ? out << "[" << P << "/" << Q << ", " << p << "/" << q << "]" :
               out << "[" << p << "/" << q << ", " << P << "/" << Q << "]";
    }
    int f;
    while ((f = P / Q) == p / q) {
        out << f;
        if (komma) (out << ".", komma = false);
        P -= f * Q;
        p -= f * q;
        P *= 10;
        p *= 10;
    }
    return out;
}

ContinuedFraction::operator int() const {
    ContinuedFraction* cf0 = copy();
    int next = cf0->has_next() ? cf0->next() : 0;
    delete cf0;
    return next;
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

Rational* Rational::copy() const {
    return new Rational(num, denom);
}

Factory::Factory(int (*next_fn)(int)) : next_fn(next_fn), n(0) {}

const bool Factory::has_next() {
    return true;
}

const int Factory::next() {
    return next_fn(n++);
}

Factory* Factory::copy() const {
    Factory* fact = new Factory(next_fn);
    fact->n = n;
    return fact;
}