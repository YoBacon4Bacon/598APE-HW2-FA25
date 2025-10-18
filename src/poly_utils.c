#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Poly create_poly(size_t degree) {
  Poly p;
  p.coeffs.resize(degree + 1, 0.0);
  return p;
}

Poly copy_poly(const Poly &p) {
  Poly copy;
  copy.coeffs = p.coeffs;
  return copy;
}

double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(const Poly &p) {
  for (int64_t i = p.coeffs.size() - 1; i >= 0; i--) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(const Poly &p, int64_t degree) {
  if (degree < 0 || (size_t)degree >= p.coeffs.size()) {
    return 0.0;
  }
  return p.coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree < 0) {
    return;
  }
  if ((size_t)degree >= p->coeffs.size()) {
    p->coeffs.resize(degree + 1, 0.0);
  }
  p->coeffs[degree] = value;
}

Poly coeff_mod(const Poly &p, double modulus) {
  Poly out;
  out.coeffs.resize(p.coeffs.size(), 0.0);
  for (size_t i = 0; i < p.coeffs.size(); i++) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      double rounded = round(p.coeffs[i]);
      double m = positive_fmod(rounded, modulus);
      out.coeffs[i] = m;
    }
  }
  return out;
}

Poly poly_add(const Poly &a, const Poly &b) {
  size_t max_size = (a.coeffs.size() > b.coeffs.size()) ? a.coeffs.size() : b.coeffs.size();
  Poly sum;
  sum.coeffs.resize(max_size, 0.0);
  for (size_t i = 0; i < a.coeffs.size(); i++) {
    sum.coeffs[i] += a.coeffs[i];
  }
  for (size_t i = 0; i < b.coeffs.size(); i++) {
    sum.coeffs[i] += b.coeffs[i];
  }
  return sum;
}

Poly poly_mul_scalar(const Poly &p, double scalar) {
  Poly res;
  res.coeffs.resize(p.coeffs.size());
  for (size_t i = 0; i < p.coeffs.size(); i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
  }
  return res;
}

Poly poly_mul(const Poly &a, const Poly &b) {
  Poly res;
  if (a.coeffs.empty() || b.coeffs.empty()) {
    res.coeffs.resize(1, 0.0);
    return res;
  }
  res.coeffs.resize(a.coeffs.size() + b.coeffs.size() - 1, 0.0);

  for (size_t i = 0; i < a.coeffs.size(); i++) {
    if (fabs(a.coeffs[i]) > 1e-9) {
      for (size_t j = 0; j < b.coeffs.size(); j++) {
        if (fabs(b.coeffs[j]) > 1e-9) {
          res.coeffs[i + j] += a.coeffs[i] * b.coeffs[j];
        }
      }
    }
  }
  return res;
}

Poly poly_rem(const Poly &num, const Poly &den) {
  // In our case `den` should always be (x^n + 1)
  size_t num_coeffs = num.coeffs.size();
  size_t n = poly_degree(den);
    Poly rem;
    rem.coeffs.resize(n, 0.0);
    size_t low_part_size = (num_coeffs < n) ? num_coeffs : n;
    for (size_t i = 0; i < low_part_size; ++i) {
        rem.coeffs[i] = num.coeffs[i];
    }
    size_t high_part_size = num_coeffs;
    for (size_t i = n; i < high_part_size; ++i) {
        rem.coeffs[i - n] -= num.coeffs[i];
    }
    return rem;
}

Poly poly_round_div_scalar(const Poly &x, double divisor) {
  Poly out;
  out.coeffs.resize(x.coeffs.size());
  assert(fabs(divisor) > 1e-9);

  for (size_t i = 0; i < x.coeffs.size(); i++) {
    double v = x.coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  return out;
}
