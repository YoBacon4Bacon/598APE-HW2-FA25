#ifndef HE_H
#define HE_H

#include "types.h"
#include <stdint.h>

typedef struct {
  PublicKey pk;
  SecretKey sk;
} KeyPair;

KeyPair keygen(size_t n, double q, const Poly &poly_mod);

Ciphertext encrypt(const PublicKey &pk, size_t n, double q, const Poly &poly_mod, double t,
                   double pt);

double decrypt(const SecretKey &sk, size_t n, double q, const Poly &poly_mod, double t,
               const Ciphertext &ct);

Poly encode_plain_integer(double t, double pt);

Ciphertext add_plain(const Ciphertext &ct, double q, double t, const Poly &poly_mod,
                     double pt);

Ciphertext add_cipher(const Ciphertext &c1, const Ciphertext &c2, double q, const Poly &poly_mod);

Ciphertext mul_plain(const Ciphertext &ct, double q, double t, const Poly &poly_mod,
                     double pt);

EvalKey evaluate_keygen(const SecretKey &sk, size_t n, double q, const Poly &poly_mod,
                        double p);

Ciphertext mul_cipher(const Ciphertext &c1, const Ciphertext &c2, double q, double t,
                      double p, const Poly &poly_mod, const EvalKey &rlk);

#endif