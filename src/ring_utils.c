#include "ring_utils.h"
#include "poly_utils.h"

void ring_add_mod(Poly *out, const Poly *x, const Poly *y, double modulus, const Poly *poly_mod) {
  Poly sum, sum_mod, quot, rem;
  poly_add(&sum, x, y);

  coeff_mod(&sum_mod, &sum, modulus);
  coeff_mod(&sum_mod, &sum_mod, modulus);

  poly_divmod(&sum_mod, poly_mod, &quot, &rem);

  coeff_mod(out, &rem, modulus);
  coeff_mod(out, out, modulus);
}

void ring_mul_mod(Poly *out, const Poly *x, const Poly *y, double modulus, const Poly *poly_mod) {
  Poly prod, prod_mod, quot, rem;
  poly_mul(&prod, x, y);

  coeff_mod(&prod_mod, &prod, modulus);
  coeff_mod(&prod_mod, &prod_mod, modulus);

  poly_divmod(&prod_mod, poly_mod, &quot, &rem);

  coeff_mod(out, &rem, modulus);
  coeff_mod(out, out, modulus);
}

void ring_mul_no_mod_q(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod) {
  Poly prod, quot;
  poly_mul(&prod, x, y);

  poly_divmod(&prod, poly_mod, &quot, out);
}

void ring_add_no_mod_q(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod) {
  Poly sum, quot;
  poly_add(&sum, x, y);

  poly_divmod(&sum, poly_mod, &quot, out);
}

void ring_mul_poly_mod(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod) {
  Poly prod, quot;
  poly_mul(&prod, x, y);
  poly_divmod(&prod, poly_mod, &quot, out);
}

void ring_add_poly_mod(Poly *out, const Poly *x, const Poly *y, const Poly *poly_mod) {
  Poly sum, quot;
  poly_add(&sum, x, y);
  poly_divmod(&sum, poly_mod, &quot, out);
}