//
// Created by Ina on 13.08.2019.
//

#ifndef CF_CONTINUEDFRACTION_H
#define CF_CONTINUEDFRACTION_H

#include <iosfwd>

class ContinuedFraction {
public:
    virtual const int next();
    virtual const bool has_next();

    // Operators:
    ContinuedFraction* operator + (ContinuedFraction &);
    ContinuedFraction* operator - (ContinuedFraction &);
    ContinuedFraction* operator / (ContinuedFraction &);
    ContinuedFraction* operator * (ContinuedFraction &);
    friend std::ostream& operator<< (std::ostream &out, ContinuedFraction &);
};

class Rational : public ContinuedFraction {
    int num, denom;
public:
    Rational(int, int);
    const int next();
    const bool has_next();
};

class Transform : public ContinuedFraction {
    virtual const bool must_feed();
    virtual void consume();
    virtual void consume(int);
};

// f(x) = (ax + b) / (cx + d)
class MoebiusTransform : public Transform {
private:
    ContinuedFraction* x;
    int a,b,c,d, term, temp;
    bool has_term = false;
public:
    MoebiusTransform(ContinuedFraction*, int, int, int, int);
    const bool must_feed();
    void consume();
    void consume(int next);
    const int next();
    const bool has_next();
};

// f(x,y) = (ax + bxy + cy + d) / (ex + fxy + gy + h)
class BivariateMoebiusTransform : public Transform {
private:
    int a,b,c,d,e,f,g,h, term, temp;
    ContinuedFraction* cfn[2];
    int cfn_index = -1;
    bool has_term = false;
    double dh, ae, cg, bf;

    int choose_cfn();
public:
    BivariateMoebiusTransform(ContinuedFraction* x, ContinuedFraction* y, int, int, int, int, int, int, int, int);
    const bool must_feed();
    void consume();
    void consume(int next);
    const int next();
    const bool has_next();
};

#endif //CF_CONTINUEDFRACTION_H
