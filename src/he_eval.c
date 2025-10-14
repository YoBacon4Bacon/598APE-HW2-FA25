#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <assert.h>
#include <math.h>
#include <algorithm>

Ciphertext add_plain(const Ciphertext &ct, double q, double t, const Poly &poly_mod,
                     double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly scaled_m = poly_mul_scalar(m, q / t);
  Poly new_c0 = ring_add_mod(ct.c0, scaled_m, q, poly_mod);

  Ciphertext result;
  result.c0 = new_c0;
  result.c1 = ct.c1;
  return result;
}

Ciphertext add_cipher(const Ciphertext &c1, const Ciphertext &c2, double q, const Poly &poly_mod) {
  Ciphertext result;
  result.c0 = ring_add_mod(c1.c0, c2.c0, q, poly_mod);
  result.c1 = ring_add_mod(c1.c1, c2.c1, q, poly_mod);
  return result;
}

Ciphertext mul_plain(const Ciphertext &ct, double q, double t, const Poly &poly_mod,
                     double pt) {
  Poly m = encode_plain_integer(t, pt);
  Ciphertext result;
  result.c0 = ring_mul_mod(ct.c0, m, q, poly_mod);
  result.c1 = ring_mul_mod(ct.c1, m, q, poly_mod);
  return result;
}

Ciphertext mul_cipher(const Ciphertext &c1, const Ciphertext &c2, double q, double t,
                      double p, const Poly &poly_mod, const EvalKey &rlk) {
  size_t n = std::max({
    c1.c0.coeffs.size(),
    c2.c0.coeffs.size(),
    c1.c1.coeffs.size(),
    c2.c1.coeffs.size()
  }) - 1;
  Poly c0_prod = ring_mul_no_mod_q(c1.c0, c2.c0, poly_mod);
  Poly c1_left = ring_mul_no_mod_q(c1.c0, c2.c1, poly_mod);
  Poly c1_right = ring_mul_no_mod_q(c1.c1, c2.c0, poly_mod);
  Poly c1_sum = ring_add_no_mod_q(c1_left, c1_right, poly_mod);
  Poly c2_prod = ring_mul_no_mod_q(c1.c1, c2.c1, poly_mod);

  Poly c0_res = create_poly(n);
  Poly c1_res = create_poly(n);
  Poly c2_res = create_poly(n);
  for (size_t i = 0; i < c0_prod.coeffs.size() && i <= n; i++) {
    if (fabs(c0_prod.coeffs[i]) > 1e-9) {
      c0_res.coeffs[i] = round(t * c0_prod.coeffs[i] / q);
    }
  }
  for (size_t i = 0; i < c1_sum.coeffs.size() && i <= n; i++) {
    if (fabs(c1_sum.coeffs[i]) > 1e-9) {
      c1_res.coeffs[i] = round(t * c1_sum.coeffs[i] / q);
    }
  }
  for (size_t i = 0; i < c2_prod.coeffs.size() && i <= n; i++) {
    if (fabs(c2_prod.coeffs[i]) > 1e-9) {
      c2_res.coeffs[i] = round(t * c2_prod.coeffs[i] / q);
    }
  }

  Poly c0_modq = coeff_mod(c0_res, q);
  Poly c1_modq = coeff_mod(c1_res, q);
  Poly c2_modq = coeff_mod(c2_res, q);

  // Relinearization
  Poly prod_b = ring_mul_no_mod_q(rlk.b, c2_modq, poly_mod);
  Poly prod_a = ring_mul_no_mod_q(rlk.a, c2_modq, poly_mod);

  Poly div_b = create_poly(n);
  Poly div_a = create_poly(n);
  for (size_t i = 0; i < prod_b.coeffs.size() && i <= n; i++) {
    double vb = prod_b.coeffs[i];
    if (fabs(vb) > 1e-9) {
      div_b.coeffs[i] = round(vb / p);
    }
  }
  for (size_t i = 0; i < prod_a.coeffs.size() && i <= n; i++) {
    double va = prod_a.coeffs[i];
    if (fabs(va) > 1e-9) {
      div_a.coeffs[i] = round(va / p);
    }
  }

  Poly c20_modq = coeff_mod(div_b, q);
  Poly c21_modq = coeff_mod(div_a, q);

  Ciphertext out;
  out.c0 = ring_add_mod(c0_modq, c20_modq, q, poly_mod);
  out.c1 = ring_add_mod(c1_modq, c21_modq, q, poly_mod);
  return out;
}