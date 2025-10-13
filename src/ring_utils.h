#ifndef RING_UTILS_H
#define RING_UTILS_H

#include "types.h"
#include <stdint.h>

void ring_add_mod(Poly *out, const Poly *x, const Poly *y, double modulus, const Poly *poly_mod);

void ring_mul_mod(Poly *out, const Poly *x, const Poly *y, double modulus, const Poly *poly_mod);

void ring_mul_no_mod_q(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod);

void ring_add_no_mod_q(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod);

void ring_add_poly_mod(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod);

void ring_mul_poly_mod(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod);

#endif