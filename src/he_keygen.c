#include "he.h"
#include "poly_random.h"
#include "poly_utils.h"
#include "ring_utils.h"

void keygen(KeyPair *out, size_t n, double q, const Poly *poly_mod) {
  SecretKey s = gen_binary_poly(n);
  Poly a = gen_uniform_poly(n, q);
  Poly e = gen_normal_poly(n, 0.0, 1.0);

  Poly neg_a, as, neg_e, b;
  poly_mul_scalar(&neg_a, &a, -1);
  ring_mul_mod(&as, &neg_a, &s, q, poly_mod);
  poly_mul_scalar(&neg_e, &e, -1);
  ring_add_mod(&b, &as, &neg_e, q, poly_mod);

  out->pk.a = a;
  out->pk.b = b;
  out->sk = s;
}

void evaluate_keygen(EvalKey *out, SecretKey sk, size_t n, double q, const Poly *poly_mod,
                     double p) {
  double new_modulus = q * p;
  Poly a = gen_uniform_poly(n, new_modulus);
  Poly e = gen_normal_poly(n, 0.0, 1.0);

  Poly s2, secret_scaled, neg_a, as, neg_e, as_nege, b_ring, b;
  poly_mul(&s2, &sk, &sk);
  poly_mul_scalar(&secret_scaled, &s2, p);

  poly_mul_scalar(&neg_a, &a, -1.0);
  ring_mul_no_mod_q(&as, &neg_a, &sk, poly_mod);
  poly_mul_scalar(&neg_e, &e, -1.0);
  ring_add_no_mod_q(&as_nege, &as, &neg_e, poly_mod);
  ring_add_no_mod_q(&b_ring, &as_nege, &secret_scaled, poly_mod);

  coeff_mod(&b, &b_ring, new_modulus);

  out->a = a;
  out->b = b;
}