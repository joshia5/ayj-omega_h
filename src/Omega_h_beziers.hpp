#ifndef OMEGA_H_BEZIERS_HPP
#define OMEGA_H_BEZIERS_HPP

#include "Omega_h_array.hpp"

namespace Omega_h {

Real intpow(const Real b, const LO e) {
  switch (e) {
  case 0: return 1.0;
  case 1: return b;
  case 2: return b*b;
  case 3: return b*b*b;
  case 4: return b*b*b*b;
  case 5: return b*b*b*b*b;
  case 6: return b*b*b*b*b*b;
  default:
    return intpow(b, e-6) * intpow(b, 6);
  }
}

OMEGA_H_DEVICE LO binomial(int n, int i);

OMEGA_H_DEVICE LO trinomial(int n, int i, int j);

OMEGA_H_DEVICE LO quadnomial(int n, int i, int j, int k);

OMEGA_H_DEVICE Real Bij(const int i, const int j, const double u,
                        const double v);

OMEGA_H_DEVICE Real Bijk(const int i, const int j, const int k, const double u,
                         const double v, const double w) {
  return intpow(u, i)*intpow(v, j)*intpow(w, k);
}

OMEGA_H_DEVICE Real Bijkl(const int i, const int j, const int k, const int l,
                          const double u, const double v, const double w,
                          const double t) {
  return intpow(u, i)*intpow(v, j)*intpow(w, k)*intpow(t, l);
}

OMEGA_H_DEVICE Real Bij(const int ij[], const double xi[]);
OMEGA_H_DEVICE Real Bijk(const int ijk[], const double xi[]);
OMEGA_H_DEVICE Real Bijkl(const int ijkl[], const double xi[]);

OMEGA_H_DEVICE static void bezierCurve(I8 P, Reals xi, Write<Real> values);
OMEGA_H_DEVICE static void bezierTriangle(I8 P, Reals xi, Write<Real> values);
OMEGA_H_DEVICE static void bezierTet(I8 P, Reals xi, Write<Real> values);

// workaround CUDA compiler bug
#ifdef OMEGA_H_USE_CUDA
__host__
#endif
    void
    assign(Mesh& a, Mesh const& b);

#define OMEGA_H_EXPL_INST_DECL(T)
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

} // namespace Omega_h

#endif
