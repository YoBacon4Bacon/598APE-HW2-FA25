#include "ring_utils.h"
#include "poly_utils.h"

Poly ring_add_mod(const Poly &x, const Poly &y, double modulus, const Poly &poly_mod) {
  Poly sum = poly_add(x, y);

  Poly sum_mod = coeff_mod(sum, modulus);

  Poly rem = poly_rem(sum_mod, poly_mod);

  Poly rem_mod = coeff_mod(rem, modulus);

  return rem_mod;
}

Poly ring_mul_mod(const Poly &x, const Poly &y, double modulus, const Poly &poly_mod) {
  Poly prod = poly_mul(x, y);

  Poly prod_mod = coeff_mod(prod, modulus);

  Poly rem = poly_rem(prod_mod, poly_mod);

  Poly rem_mod = coeff_mod(rem, modulus);

  return rem_mod;
}

Poly ring_mul_no_mod_q(const Poly &x, const Poly &y, const Poly &poly_mod) {
  Poly prod = poly_mul(x, y);

  Poly rem = poly_rem(prod, poly_mod);

  return rem;
}

Poly ring_add_no_mod_q(const Poly &x, const Poly &y, const Poly &poly_mod) {
  Poly sum = poly_add(x, y);

  Poly rem = poly_rem(sum, poly_mod);

  return rem;
}

Poly ring_mul_poly_mod(const Poly &x, const Poly &y, const Poly &poly_mod) {
  Poly prod = poly_mul(x, y);
  Poly rem = poly_rem(prod, poly_mod);
  return rem;
}

Poly ring_add_poly_mod(const Poly &x, const Poly &y, const Poly &poly_mod) {
  Poly sum = poly_add(x, y);
  Poly rem = poly_rem(sum, poly_mod);
  return rem;
}