#include "he.h"
#include "poly_random.h"
#include "poly_utils.h"
#include "ring_utils.h"

Poly encode_plain_integer(double t, double pt) {
  Poly m = create_poly();
  double v = positive_fmod(pt, t);
  m.coeffs[0] = v;
  return m;
}

void encrypt(Ciphertext *out, PublicKey pk, size_t n, double q, const Poly *poly_mod, double t,
             double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly scaled_m;
  poly_mul_scalar(&scaled_m, &m, floor(q / t));
  Poly e1 = gen_normal_poly(n, 0.0, 1.0);
  Poly e2 = gen_normal_poly(n, 0.0, 1.0);
  Poly u = gen_binary_poly(n);

  Poly bu, bu_e1, au;
  ring_mul_mod(&bu, &pk.b, &u, q, poly_mod);
  ring_add_mod(&bu_e1, &bu, &e1, q, poly_mod);
  ring_add_mod(&out->c0, &bu_e1, &scaled_m, q, poly_mod);

  ring_mul_mod(&au, &pk.a, &u, q, poly_mod);
  ring_add_mod(&out->c1, &au, &e2, q, poly_mod);
}