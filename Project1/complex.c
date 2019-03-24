#include "complex.h"

#include <math.h>

Complex cadd(Complex lhs, Complex rhs) {
    Complex res = {lhs.re + rhs.re, lhs.im + rhs.im};
    return res;
}

Complex csub(Complex lhs, Complex rhs) {
    Complex res = {lhs.re - rhs.re, lhs.im - rhs.im};
    return res;
}

Complex cmul(Complex lhs, Complex rhs) {
    Complex res =  {lhs.re*rhs.re - lhs.im*rhs.im,
                    lhs.re*rhs.im + lhs.im*rhs.re};
    return res;
}

Complex cdiv(Complex lhs, Complex rhs) {
    double rhs2 = rhs.re*rhs.re + rhs.im*rhs.im;
    Complex res = {
        (lhs.re*rhs.re + lhs.im*rhs.im) / rhs2,
        (lhs.re*rhs.im - lhs.im*rhs.re) / rhs2
    };
    return res;
}

Complex cconj(Complex z) {
    Complex res = {z.re, -z.im};
    return res;
}

double cabs_(Complex z) {
    return sqrt(z.re*z.re + z.im*z.im);
}