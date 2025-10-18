#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <math.h>
#include <stdlib.h>

double decrypt(const SecretKey &sk, size_t n, double q, const Poly &poly_mod, double t,
               const Ciphertext &ct) {
  Poly c1s = ring_mul_mod(ct.c1, sk, q, poly_mod);
  Poly scaled_pt = ring_add_mod(c1s, ct.c0, q, poly_mod);

  if (!scaled_pt.coeffs.empty() && fabs(scaled_pt.coeffs[0]) > 1e-9) {
    double v = round(scaled_pt.coeffs[0]);
    double result = round(t * v / q);
    return positive_fmod(result, t);
  }
  
  return 0.0;
}
