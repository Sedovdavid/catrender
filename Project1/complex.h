#ifndef COMPLEX_H
#define COMPLEX_H

typedef struct {
    double re;
    double im;
} Complex;

Complex cadd(Complex, Complex);
Complex csub(Complex, Complex);
Complex cmul(Complex, Complex);
Complex cdiv(Complex, Complex);

Complex cconj(Complex);
double cabs_(Complex);
#endif