#ifndef __GPPCOMPLEX
#define __GPPCOMPLEX

#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <complex>
#include <omp.h>
#include <ctime>
//#include <chrono>
#include <stdio.h>

class GPPComplex {
    private : 
    double re;
    double im;

public:
explicit GPPComplex () {
    re = 0.00;
    im = 0.00;
}


explicit GPPComplex(const double& x, const double& y) {
    re = x;
    im = y;
}

GPPComplex(const GPPComplex& src) {
    re = src.re;
    im = src.im;
}

GPPComplex& operator =(const GPPComplex& src) {
    re = src.re;
    im = src.im;

    return *this;
}

GPPComplex& operator +(const GPPComplex& src) {
    re += src.re;
    im += src.im;

    return *this;
}

GPPComplex& operator +=(const GPPComplex& src) {
    re = src.re + this->re;
    im = src.im + this->im;

    return *this;
}

GPPComplex& operator -=(const GPPComplex& src) {
    re = src.re - this->re;
    im = src.im - this->im;

    return *this;
}

GPPComplex& operator -() {
    re = -this->re;
    im = -this->im;

    return *this;
}

GPPComplex& operator ~() {
    return *this;
}

void print() const {
    printf("( %f, %f) ", this->re, this->im);
    printf("\n");
}

double abs(const GPPComplex& src) {

    double re_this = src.re * src.re;
    double im_this = src.im * src.im;

    double result = sqrt(re_this+im_this);

    return result;

}

double get_real() const
{
    return this->re;
}

double get_imag() const
{
    return this->im;
}

void set_real(double val)
{
    this->re = val;
}

void set_imag(double val) 
{
    this->im = val;
}


    friend inline GPPComplex GPPComplex_square(GPPComplex& src) ;
    friend inline GPPComplex GPPComplex_conj(const GPPComplex& src) ;
    friend inline GPPComplex GPPComplex_product(const GPPComplex& a, const GPPComplex& b) ;
    friend inline double GPPComplex_abs(const GPPComplex& src) ;
    friend inline GPPComplex GPPComplex_mult(GPPComplex& a, double b, double c) ;
    friend inline GPPComplex GPPComplex_mult(const GPPComplex& a, double b) ;
    friend inline void GPPComplex_fma(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) ;
    friend inline void GPPComplex_fms(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) ;
    friend inline GPPComplex doubleMinusGPPComplex(const double &a, GPPComplex& src) ;
    friend inline GPPComplex doublePlusGPPComplex(double a, GPPComplex& src) ;
    friend inline double GPPComplex_real( const GPPComplex& src) ;
    friend inline double GPPComplex_imag( const GPPComplex& src) ;
};

    inline GPPComplex GPPComplex_square(GPPComplex& src) ;
    inline GPPComplex GPPComplex_conj(const GPPComplex& src) ;
    inline GPPComplex GPPComplex_product(const GPPComplex& a, const GPPComplex& b) ;
    inline double GPPComplex_abs(const GPPComplex& src) ;
    inline GPPComplex GPPComplex_mult(GPPComplex& a, double b, double c) ;
    inline GPPComplex GPPComplex_mult(const GPPComplex& a, double b) ;
    inline void GPPComplex_fma(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) ;
    inline void GPPComplex_fms(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) ;

//Inline functions have to be defined in the same file as the declaration

/*
 * Return the square of a complex number 
 */
GPPComplex GPPComplex_square(GPPComplex& src) {
    double re_this = src.re ;
    double im_this = src.im ;

    GPPComplex result(re_this*re_this - im_this*im_this, 2*re_this*im_this);

    return result;
}

/*
 * Return the conjugate of a complex number 
 */
GPPComplex GPPComplex_conj(const GPPComplex& src) {

double re_this = src.re;
double im_this = -1 * src.im;

GPPComplex result(re_this, im_this);
return result;

}


/*
 * Return the product of 2 complex numbers 
 */
inline GPPComplex GPPComplex_product(const GPPComplex& a, const GPPComplex& b) {

    double re_this = a.re * b.re - a.im*b.im ;
    double im_this = a.re * b.im + a.im*b.re ;

    GPPComplex result(re_this, im_this);
    return result;
}

/*
 * Return the absolute of a complex number 
 */
double GPPComplex_abs(const GPPComplex& src) {
    double re_this = src.re * src.re;
    double im_this = src.im * src.im;

    double result = sqrt(re_this+im_this);
    return result;
}

/*
 *  result = a * b * c (a = complex ; b,c = double) 
 */
GPPComplex GPPComplex_mult(GPPComplex& a, double b, double c) {

    GPPComplex result(a.re * b * c, a.im * b * c);
    return result;

}

/*
 * Return the complex number c = a * b (a is complex, b is double) 
 */
GPPComplex GPPComplex_mult(const GPPComplex& a, double b) {

   GPPComplex result(a.re*b, a.im*b);
   return result;

}

/*
 * Return the complex number a += b * c  
 */
void GPPComplex_fma(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) {
    double re_this = b.re * c.re - b.im*c.im ;
    double im_this = b.re * c.im + b.im*c.re ;

    GPPComplex mult_result(re_this, im_this);

    a.re += mult_result.re;
    a.im += mult_result.im;
}

/*
 * Return the complex number a -= b * c  
 */
void GPPComplex_fms(GPPComplex& a, const GPPComplex& b, const GPPComplex& c) {
    double re_this = b.re * c.re - b.im*c.im ;
    double im_this = b.re * c.im + b.im*c.re ;

    GPPComplex mult_result(re_this, im_this);

    a.re -= mult_result.re;
    a.im -= mult_result.im;
}


GPPComplex doubleMinusGPPComplex(const double &a, GPPComplex& src) {
    GPPComplex result(a - src.re, 0 - src.im);
    return result;
}

GPPComplex doublePlusGPPComplex(double a, GPPComplex& src) {
    GPPComplex result(a + src.re, 0 + src.im);
    return result;
}

double GPPComplex_real( const GPPComplex& src) {
    return src.re;
}

double GPPComplex_imag( const GPPComplex& src) {
    return src.im;
}

#endif
