#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "types.h"
#include <math.h>
#include <stdint.h>

double positive_fmod(double x, double m);

int64_t poly_degree(Poly p);

double get_coeff(Poly p, int64_t degree);

void set_coeff(Poly *p, int64_t degree, double value);

void coeff_mod(Poly *out, const Poly *p, double modulus);

void poly_add(Poly *out, const Poly *a, const Poly *b);

void poly_mul_scalar(Poly *out, const Poly *p, double scalar);

void poly_mul(Poly *out, const Poly *a, const Poly *b);

void poly_divmod(const Poly *numerator, const Poly *denominator, Poly *quotient,
                 Poly *remainder);

Poly create_poly(void);

#endif