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
  assert(poly_degree(den) > 0 || fabs(get_coeff(den, 0)) > 1e-9);

  size_t ndeg = poly_degree(num);
  size_t ddeg = poly_degree(den);

  Poly rem = copy_poly(num);

  if (ndeg < ddeg) {
    return rem;
  }

  double d_lead = get_coeff(den, ddeg);
  assert(fabs(d_lead) > 1e-9);

  for (int64_t k = ndeg - ddeg; k >= 0; --k) {
    int64_t target_deg = ddeg + k;
    double r_coeff = get_coeff(rem, target_deg);
    double coeff = trunc(round(r_coeff) / round(d_lead));

    for (size_t i = 0; i < den.coeffs.size(); i++) {
      if (fabs(den.coeffs[i]) > 1e-9) {
        int64_t deg = i + k;
        rem.coeffs[deg] -= coeff * den.coeffs[i];
      }
    }
  }

  assert(poly_degree(rem) < poly_degree(den));

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
