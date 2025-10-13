#ifndef HE_H
#define HE_H

#include "types.h"
#include <stdint.h>

typedef struct {
  PublicKey pk;
  SecretKey sk;
} KeyPair;

void keygen(KeyPair *out, size_t n, double q, const Poly *poly_mod);

void encrypt(Ciphertext *out, PublicKey pk, size_t n, double q, const Poly *poly_mod, double t,
             double pt);

double decrypt(SecretKey sk, size_t n, double q, const Poly *poly_mod, double t,
               Ciphertext ct);

Poly encode_plain_integer(double t, double pt);

void add_plain(Ciphertext *out, Ciphertext ct, double q, double t, const Poly *poly_mod,
               double pt);

void add_cipher(Ciphertext *out, Ciphertext c1, Ciphertext c2, double q, const Poly *poly_mod);

void mul_plain(Ciphertext *out, Ciphertext ct, double q, double t, const Poly *poly_mod,
               double pt);

void evaluate_keygen(EvalKey *out, SecretKey sk, size_t n, double q, const Poly *poly_mod,
                     double p);

void mul_cipher(Ciphertext *out, Ciphertext c1, Ciphertext c2, double q, double t,
                double p, const Poly *poly_mod, EvalKey rlk);

#endif