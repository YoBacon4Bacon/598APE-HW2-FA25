#ifndef TYPES_H
#define TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <vector>

typedef struct {
  std::vector<double> coeffs;  // resizable array
} Poly;

typedef struct {
  Poly b;
  Poly a;
} PublicKey;

typedef Poly SecretKey;

typedef struct {
  Poly c0;
  Poly c1;
} Ciphertext;

typedef struct {
  Poly c0;
  Poly c1;
  Poly c2;
} Ciphertext3;

typedef struct {
  Poly a;
  Poly b;
} EvalKey;

#endif