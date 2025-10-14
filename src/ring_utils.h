#ifndef RING_UTILS_H
#define RING_UTILS_H

#include "types.h"
#include <stdint.h>

Poly ring_add_mod(const Poly &x, const Poly &y, double modulus, const Poly &poly_mod);

Poly ring_mul_mod(const Poly &x, const Poly &y, double modulus, const Poly &poly_mod);

Poly ring_mul_no_mod_q(const Poly &x, const Poly &y, const Poly &poly_mod);

Poly ring_add_no_mod_q(const Poly &x, const Poly &y, const Poly &poly_mod);

Poly ring_add_poly_mod(const Poly &x, const Poly &y, const Poly &poly_mod);

Poly ring_mul_poly_mod(const Poly &x, const Poly &y, const Poly &poly_mod);

#endif