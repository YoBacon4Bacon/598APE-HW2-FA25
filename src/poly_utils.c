#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Poly create_poly(size_t degree) {
  assert(degree < MAX_POLY_DEGREE);
  Poly p;
  for (int i = 0; i <= degree; i++) {
    p.coeffs[i] = 0.0;
  }
  p.max_degree = degree;
  return p;
}

Poly copy_poly(const Poly &p) {
  Poly copy;
  for (size_t i = 0; i <= p.max_degree; i++) {
    copy.coeffs[i] = p.coeffs[i];
  }
  copy.max_degree = p.max_degree;
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
  for (int64_t i = p.max_degree; i >= 0; i--) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(const Poly &p, int64_t degree) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p.coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  p->coeffs[degree] = value;
  p->max_degree = degree > p->max_degree ? degree : p->max_degree;
}

Poly coeff_mod(const Poly &p, double modulus) {
  Poly out = create_poly(p.max_degree);
  for (size_t i = 0; i <= p.max_degree; i++) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      double rounded = round(p.coeffs[i]);
      double m = positive_fmod(rounded, modulus);
      out.coeffs[i] = m;
    }
  }
  return out;
}

Poly poly_add(const Poly &a, const Poly &b) {
  size_t max_deg = (a.max_degree > b.max_degree) ? a.max_degree : b.max_degree;
  Poly sum = create_poly(max_deg);
  for (size_t i = 0; i <= a.max_degree; i++) {
    sum.coeffs[i] += a.coeffs[i];
  }
  for (size_t i = 0; i <= b.max_degree; i++) {
    sum.coeffs[i] += b.coeffs[i];
  }
  return sum;
}

Poly poly_mul_scalar(const Poly &p, double scalar) {
  Poly res = create_poly(p.max_degree);
  for (size_t i = 0; i <= p.max_degree; i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
  }
  return res;
}

Poly poly_mul(const Poly &a, const Poly &b) {
  Poly res = create_poly(a.max_degree + b.max_degree);

  for (size_t i = 0; i <= a.max_degree; i++) {
    if (fabs(a.coeffs[i]) > 1e-9) {
      for (size_t j = 0; j <= b.max_degree; j++) {
        if (fabs(b.coeffs[j]) > 1e-9) {
          assert(i + j < MAX_POLY_DEGREE);
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

    for (size_t i = 0; i <= den.max_degree; i++) {
      if (fabs(den.coeffs[i]) > 1e-9) {
        int64_t deg = i + k;
        assert(deg < MAX_POLY_DEGREE);
        rem.coeffs[deg] -= coeff * den.coeffs[i];
      }
    }
  }

  assert(poly_degree(rem) < poly_degree(den));

  // why not
  rem.max_degree = poly_degree(rem);

  return rem;
}

Poly poly_round_div_scalar(const Poly &x, double divisor) {
  Poly out = create_poly(x.max_degree);
  assert(fabs(divisor) > 1e-9);

  for (size_t i = 0; i <= x.max_degree; i++) {
    double v = x.coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  return out;
}
