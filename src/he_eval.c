#include "he.h"
#include "poly_utils.h"
#include "ring_utils.h"
#include <assert.h>
#include <math.h>

void add_plain(Ciphertext *out, Ciphertext ct, double q, double t, const Poly *poly_mod,
               double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly scaled_m;
  poly_mul_scalar(&scaled_m, &m, q / t);
  ring_add_mod(&out->c0, &ct.c0, &scaled_m, q, poly_mod);
  out->c1 = ct.c1;
}

void add_cipher(Ciphertext *out, Ciphertext c1, Ciphertext c2, double q, const Poly *poly_mod) {
  ring_add_mod(&out->c0, &c1.c0, &c2.c0, q, poly_mod);
  ring_add_mod(&out->c1, &c1.c1, &c2.c1, q, poly_mod);
}

void mul_plain(Ciphertext *out, Ciphertext ct, double q, double t, const Poly *poly_mod,
               double pt) {
  Poly m = encode_plain_integer(t, pt);
  ring_mul_mod(&out->c0, &ct.c0, &m, q, poly_mod);
  ring_mul_mod(&out->c1, &ct.c1, &m, q, poly_mod);
}

void mul_cipher(Ciphertext *out, Ciphertext c1, Ciphertext c2, double q, double t,
                double p, const Poly *poly_mod, EvalKey rlk) {
  Poly c0_prod, c1_left, c1_right, c1_sum, c2_prod;
  ring_mul_no_mod_q(&c0_prod, &c1.c0, &c2.c0, poly_mod);
  ring_mul_no_mod_q(&c1_left, &c1.c0, &c2.c1, poly_mod);
  ring_mul_no_mod_q(&c1_right, &c1.c1, &c2.c0, poly_mod);
  ring_add_no_mod_q(&c1_sum, &c1_left, &c1_right, poly_mod);
  ring_mul_no_mod_q(&c2_prod, &c1.c1, &c2.c1, poly_mod);

  Poly c0_res = create_poly();
  Poly c1_res = create_poly();
  Poly c2_res = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    if (fabs(c0_prod.coeffs[i]) > 1e-9) {
      c0_res.coeffs[i] = round(t * c0_prod.coeffs[i] / q);
    }
    if (fabs(c1_sum.coeffs[i]) > 1e-9) {
      c1_res.coeffs[i] = round(t * c1_sum.coeffs[i] / q);
    }
    if (fabs(c2_prod.coeffs[i]) > 1e-9) {
      c2_res.coeffs[i] = round(t * c2_prod.coeffs[i] / q);
    }
  }

  Poly c0_modq, c1_modq, c2_modq;
  coeff_mod(&c0_modq, &c0_res, q);
  coeff_mod(&c1_modq, &c1_res, q);
  coeff_mod(&c2_modq, &c2_res, q);

  // Relinearization
  Poly prod_b, prod_a;
  ring_mul_no_mod_q(&prod_b, &rlk.b, &c2_modq, poly_mod);
  ring_mul_no_mod_q(&prod_a, &rlk.a, &c2_modq, poly_mod);

  Poly div_b = create_poly();
  Poly div_a = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double vb = prod_b.coeffs[i];
    double va = prod_a.coeffs[i];
    if (fabs(vb) > 1e-9) {
      div_b.coeffs[i] = round(vb / p);
    }
    if (fabs(va) > 1e-9) {
      div_a.coeffs[i] = round(va / p);
    }
  }

  Poly c20_modq, c21_modq;
  coeff_mod(&c20_modq, &div_b, q);
  coeff_mod(&c21_modq, &div_a, q);

  ring_add_mod(&out->c0, &c0_modq, &c20_modq, q, poly_mod);
  ring_add_mod(&out->c1, &c1_modq, &c21_modq, q, poly_mod);
}