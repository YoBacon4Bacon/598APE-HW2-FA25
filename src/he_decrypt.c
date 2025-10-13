#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <math.h>
#include <stdlib.h>

double decrypt(SecretKey sk, size_t n, double q, const Poly *poly_mod, double t,
                Ciphertext ct) {
  Poly c1s, scaled_pt;
  ring_mul_mod(&c1s, &ct.c1, &sk, q, poly_mod);
  ring_add_mod(&scaled_pt, &c1s, &ct.c0, q, poly_mod);

  Poly dec = create_poly();

  for (int64_t i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(scaled_pt.coeffs[i]) > 1e-9) {
      double v = round(scaled_pt.coeffs[i]);
      double result = round(t * v / q);
      dec.coeffs[i] = positive_fmod(result, t);
    }
  }

  return round(get_coeff(dec, 0));
}
