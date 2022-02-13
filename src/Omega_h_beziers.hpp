#ifndef OMEGA_H_BEZIERS_HPP
#define OMEGA_H_BEZIERS_HPP

#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_few.hpp"

namespace Omega_h {

constexpr OMEGA_H_INLINE Real const_factorial(Int N) {
  switch (N) {
    case 0:
      return 1.0;
    case 1:
      return 1.0;
    case 2:
      return 2.0;
    case 3:
      return 6.0;
    case 4:
      return 24.0;
    case 5:
      return 120.0;
    case 6:
      return 720.0;
    case 7:
      return 5040.0;
    case 8:
      return 40320.0;
    case 9:
      return 362880.0;
    case 10:
      return 3628800.0;
  }
  return -1.0;
}

OMEGA_H_INLINE Real B0_quad(Real u) {
  return (1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B1_quad(Real u) {
  return 2.0*u*(1.0-u);
}

OMEGA_H_INLINE Real B2_quad(Real u) {
  return u*u;
}

OMEGA_H_INLINE Real B0_cube(Real u) {
  return (1.0-u)*(1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B1_cube(Real u) {
  return 3.0*u*(1.0-u)*(1.0-u);
}

OMEGA_H_INLINE Real B2_cube(Real u) {
  return 3.0*(1.0-u)*u*u;
}

OMEGA_H_INLINE Real B3_cube(Real u) {
  return u*u*u;
}

Real B0_quart(Real u);
Real B1_quart(Real u);
Real B2_quart(Real u);
Real B3_quart(Real u);
Real B4_quart(Real u);

Real B00_quad(Real u, Real v);
Real B10_quad(Real u, Real v);
Real B20_quad(Real u, Real v);
Real B11_quad(Real u, Real v);
Real B02_quad(Real u, Real v);
Real B01_quad(Real u, Real v);

Real B00_cube(Real u, Real v);
Real B10_cube(Real u, Real v);
Real B20_cube(Real u, Real v);
Real B30_cube(Real u, Real v);
Real B21_cube(Real u, Real v);
Real B12_cube(Real u, Real v);
Real B03_cube(Real u, Real v);
Real B02_cube(Real u, Real v);
Real B01_cube(Real u, Real v);
Real B11_cube(Real u, Real v);

Real B00_quart(Real u, Real v);
Real B10_quart(Real u, Real v);
Real B20_quart(Real u, Real v);
Real B30_quart(Real u, Real v);
Real B40_quart(Real u, Real v);
Real B31_quart(Real u, Real v);
Real B22_quart(Real u, Real v);
Real B13_quart(Real u, Real v);
Real B04_quart(Real u, Real v);
Real B03_quart(Real u, Real v);
Real B02_quart(Real u, Real v);
Real B01_quart(Real u, Real v);
Real B11_quart(Real u, Real v);
Real B21_quart(Real u, Real v);
Real B12_quart(Real u, Real v);

// Babushka chen points
constexpr OMEGA_H_INLINE Real xi_1_quad() {
  return 0.5;
}

constexpr OMEGA_H_INLINE Real xi_1_cube() {
  return 0.2748043;
}
constexpr OMEGA_H_INLINE Real xi_2_cube() {
  return 0.7251957;
}

constexpr OMEGA_H_INLINE Real xi_1_quart() {
  return 0.1693976;
}
constexpr OMEGA_H_INLINE Real xi_2_quart() {
  return 0.5;
}
constexpr OMEGA_H_INLINE Real xi_3_quart() {
  return 0.8306024;
}

constexpr OMEGA_H_INLINE Real xi_1_quint() {
  return 0.1133573;
}
constexpr OMEGA_H_INLINE Real xi_2_quint() {
  return 0.3568239;
}
constexpr OMEGA_H_INLINE Real xi_3_quint() {
  return 0.6431761;
}
constexpr OMEGA_H_INLINE Real xi_4_quint() {
  return 0.8866427;
}

constexpr OMEGA_H_INLINE Real xi_1_hex() {
  return 0.0805979;
}
constexpr OMEGA_H_INLINE Real xi_2_hex() {
  return 0.2650895;
}
constexpr OMEGA_H_INLINE Real xi_3_hex() {
  return 0.5;
}
constexpr OMEGA_H_INLINE Real xi_4_hex() {
  return 0.7349105;
}
constexpr OMEGA_H_INLINE Real xi_5_hex() {
  return 0.9194021;
}

OMEGA_H_INLINE Vector<2> xi_11_cube() {
  return vector_2(1.0/3.0, 1.0/3.0);
}
OMEGA_H_INLINE Vector<2> xi_11_quart() {
  return vector_2(0.22088805, 0.22088805);
}
OMEGA_H_INLINE Vector<2> xi_21_quart() {
  return vector_2(0.5582239, 0.22088805);
}
OMEGA_H_INLINE Vector<2> xi_12_quart() {
  return vector_2(0.22088805, 0.5582239);
}
OMEGA_H_INLINE Vector<3> xi_111_quart() {
  return vector_3(1.0/4.0, 1.0/4.0, 1.0/4.0);
}

/*
//TODO limititations for generalized order : 1. these fns will need to be
//dynamically sized, 2. lin eqn form and solve on gpu 
Reals curve_bezier_pts(LO const P) {
//OMEGA_H_INLINE Reals curve_bezier_pts(LO const P) {
  switch (P) {
    case 2: {
      return Reals({xi_1_quad()});
    }
    case 3: {
      return Reals({xi_1_cube(),xi_2_cube()});
    }
    case 4: {
      return Reals({xi_1_quart(),xi_2_quart(),xi_3_quart()});
    }
    case 5: {
      return Reals({xi_1_quint(),xi_2_quint(),xi_3_quint(),xi_4_quint()});
    }
    case 6: {
      return Reals({xi_1_hex(),xi_2_hex(),xi_3_hex(),xi_4_hex(),xi_5_hex()});
    }
  }
  OMEGA_H_NORETURN(Reals());
}

Reals triangle_bezier_pts(LO const P) {
//OMEGA_H_INLINE Reals triangle_bezier_pts(LO const P) {
  switch (P) {
    case 3: {
      return xi_11_cube();
    }
    case 4: {
      return Read<Real>(concat<Real>(concat<Real>(xi_11_quart(), xi_21_quart()), xi_12_quart()));
    }
    case 5: {
      return Reals({});
    }
    case 6: {
      return Reals({});
    }
  }
  OMEGA_H_NORETURN(Reals());
}

OMEGA_H_INLINE Reals tet_bezier_pts(LO const P) {
  switch (P) {
    case 4: {
      return xi_111_quart();
    }
    case 5: {
      return Reals({});
    }
    case 6: {
      return Reals({});
    }
  }
  OMEGA_H_NORETURN(Reals());
}
*/

OMEGA_H_INLINE Real
Bijk(LO const P, LO const i, LO const j, LO const k, Real const u,
     Real const v, Real const w) noexcept {
  LO const l = P - i - j - k;
  OMEGA_H_CHECK(l >= 0);
  Real const t = 1.0 - u - v - w;
  OMEGA_H_CHECK((t >= 0.0) && (t <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*const_factorial(P);
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*std::pow(w,k);
  resultant = resultant*std::pow(t,l);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);
  resultant = resultant/const_factorial(k);
  resultant = resultant/const_factorial(l);

  return resultant;
}

OMEGA_H_INLINE Real 
Bij(LO const P, LO const i, LO const j, Real const u, Real const v) noexcept {
  LO const k = P - i - j;
  OMEGA_H_CHECK((k >= 0) && (k <= P));
  Real const w = 1.0 - u - v;
  OMEGA_H_CHECK((w >= 0.0) && (w <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*const_factorial(P);
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*std::pow(w,k);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);
  resultant = resultant/const_factorial(k);

  return resultant;
}

OMEGA_H_INLINE Real 
Bi(LO const P, LO const i, Real const u) noexcept {
  LO const j = P - i;
  OMEGA_H_CHECK((j >= 0) && (j <= P));
  Real const v = 1.0 - u;
  OMEGA_H_CHECK((v >= 0.0) && (v <= 1.0));
  Real resultant = 1.0;
  resultant = resultant*std::pow(u,i);
  resultant = resultant*std::pow(v,j);
  resultant = resultant*const_factorial(P);
  resultant = resultant/const_factorial(i);
  resultant = resultant/const_factorial(j);

  return resultant;
}

void elevate_curve_order_2to3(Mesh* mesh);
void elevate_curve_order_3to4(Mesh* mesh);
void elevate_curve_order_4to5(Mesh* mesh);
void elevate_curve_order_5to6(Mesh* mesh);

void calc_quad_ctrlPts_from_interpPts(Mesh *mesh);

OMEGA_H_INLINE Few<Real,4> cubic_noKeyEdge_xi_values(
    LO old_vert, LO v0, LO v1, LO v2, LO old_edge, LO e0, LO e1, LO e2) {

  Few<Real, 4> xi1_xi2;
  if (old_vert == v0) {
    OMEGA_H_CHECK(old_edge == e1);
    xi1_xi2[0] = 0.13740215;
    xi1_xi2[1] = 0.13740215;
    xi1_xi2[2] = 0.362597849;
    xi1_xi2[3] = 0.362597849;
  }
  else if (old_vert == v1) {
    OMEGA_H_CHECK(old_edge == e2);
    xi1_xi2[0] = 0.7251957;
    xi1_xi2[1] = 0.13740215;
    xi1_xi2[2] = 0.2748043;
    xi1_xi2[3] = 0.362597849;
  }
  else if (old_vert == v2) {
    OMEGA_H_CHECK(old_edge == e0);
    xi1_xi2[0] = 0.13740215;
    xi1_xi2[1] = 0.7251957;
    xi1_xi2[2] = 0.362597849;
    xi1_xi2[3] = 0.2748043;
  }
  else {
  }

  return xi1_xi2;
}

OMEGA_H_INLINE Few<Real,4> cubic_faceSplittingEdge_xi_values
  (LO const old_noKey_vert, LO const old_v0, LO const old_v1, LO const old_v2,
   LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const new_v0_f0, LO const new_v1_f0, LO const new_v2_f0,
   LO const new_v0_f1, LO const new_v1_f1, LO const new_v2_f1,
   LO const oldf_v0_new, LO const oldf_v1_new, LO const oldf_v2_new) {

  Few<Real, 4> xi1_xi2;
  I8 should_swap = -1;
  if (old_noKey_vert == old_v0) {
    OMEGA_H_CHECK(old_key_edge == old_e1);
    xi1_xi2[0] = 1.0/2.0;
    xi1_xi2[1] = 1.0/6.0;
    xi1_xi2[2] = 1.0/6.0;
    xi1_xi2[3] = 1.0/2.0;

    if ((oldf_v1_new == new_v0_f0) || (oldf_v1_new == new_v1_f0) || 
        (oldf_v1_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v1_new == new_v0_f1) || (oldf_v1_new == new_v1_f1) ||
                     (oldf_v1_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else if (old_noKey_vert == old_v1) {
    OMEGA_H_CHECK(old_key_edge == old_e2);
    xi1_xi2[0] = 1.0/3.0;
    xi1_xi2[1] = 1.0/2.0;
    xi1_xi2[2] = 1.0/3.0;
    xi1_xi2[3] = 1.0/6.0;
    if ((oldf_v2_new == new_v0_f0) || (oldf_v2_new == new_v1_f0) ||
        (oldf_v2_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v2_new == new_v0_f1) || (oldf_v2_new == new_v1_f1) ||
                     (oldf_v2_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else if (old_noKey_vert == old_v2) {
    OMEGA_H_CHECK(old_key_edge == old_e0);
    xi1_xi2[0] = 1.0/6.0;
    xi1_xi2[1] = 1.0/3.0;
    xi1_xi2[2] = 1.0/2.0;
    xi1_xi2[3] = 1.0/3.0;
    if ((oldf_v0_new == new_v0_f0) || (oldf_v0_new == new_v1_f0) ||
        (oldf_v0_new == new_v2_f0)) {
      should_swap = -1;
    }
    else {
      OMEGA_H_CHECK ((oldf_v0_new == new_v0_f1) || (oldf_v0_new == new_v1_f1) ||
                     (oldf_v0_new == new_v2_f1));
      should_swap = 1;
    }
  }
  else {
  }
  if (should_swap == 1) {
    swap2(xi1_xi2[0], xi1_xi2[2]);
    swap2(xi1_xi2[1], xi1_xi2[3]);
  }

  return xi1_xi2;
}

OMEGA_H_INLINE Vector<3> cubic_region_xi_values
  (LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const old_e3, LO const old_e4, LO const old_e5) {
  Vector<3> p11;
  if (old_key_edge == old_e0) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e1) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e2) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/3.0;
  }
  else if (old_key_edge == old_e3) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/6.0;
  }
  else if (old_key_edge == old_e4) {
    p11[0] = 1.0/6.0;
    p11[1] = 1.0/3.0;
    p11[2] = 1.0/6.0;
  }
  else if (old_key_edge == old_e5) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/6.0;
  }
  else {
  }

  return p11;
}

OMEGA_H_INLINE LO edge_is_flip(LO const e0v0, LO const e0v1, LO const v0,
                               LO const v1) {
  LO is_flip = -1;
  if ((e0v0 == v1) && (e0v1 == v0)) {
    is_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
  }
  return is_flip;
}

OMEGA_H_INLINE Int n_internal_ctrlPts(Int edim, Int max_order) {
  OMEGA_H_CHECK(max_order > 0);
  OMEGA_H_CHECK(edim >= 0);
  if (edim == 0) {
    return 1;
  }
  else if (edim == 1) {
    return max_order-1;
  }
  else if (edim == 2) {
    return ((max_order-1)*(max_order-2))/2;
  }
  else if (edim == 3) {
    return ((max_order-1)*(max_order-2)*(max_order-3))/6;
  }
  else {
    return -1;
  }
  return -1;
}

OMEGA_H_DEVICE Vector<3> rgn_parametricToParent_3d(
    LO const order, LO const old_rgn, LOs old_ev2v, LOs old_rv2v,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Vector<3> nodePt, LOs old_re2e, LOs old_rf2f) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  auto old_rgn_v0 = old_rv2v[old_rgn*4 + 0];
  auto old_rgn_v1 = old_rv2v[old_rgn*4 + 1];
  auto old_rgn_v2 = old_rv2v[old_rgn*4 + 2];
  auto old_rgn_v3 = old_rv2v[old_rgn*4 + 3];
  auto c000 = get_vector<3>(old_vertCtrlPts, old_rgn_v0);
  auto c300 = get_vector<3>(old_vertCtrlPts, old_rgn_v1);
  auto c030 = get_vector<3>(old_vertCtrlPts, old_rgn_v2);
  auto c003 = get_vector<3>(old_vertCtrlPts, old_rgn_v3);

  auto old_rgn_e0 = old_re2e[old_rgn*6 + 0];
  auto old_rgn_e1 = old_re2e[old_rgn*6 + 1];
  auto old_rgn_e2 = old_re2e[old_rgn*6 + 2];
  auto old_rgn_e3 = old_re2e[old_rgn*6 + 3];
  auto old_rgn_e4 = old_re2e[old_rgn*6 + 4];
  auto old_rgn_e5 = old_re2e[old_rgn*6 + 5];

  auto old_rgn_f0 = old_rf2f[old_rgn*4 + 0];
  auto old_rgn_f1 = old_rf2f[old_rgn*4 + 1];
  auto old_rgn_f2 = old_rf2f[old_rgn*4 + 2];
  auto old_rgn_f3 = old_rf2f[old_rgn*4 + 3];
  auto c110 = get_vector<3>(old_faceCtrlPts, old_rgn_f0);
  auto c101 = get_vector<3>(old_faceCtrlPts, old_rgn_f1);
  auto c111 = get_vector<3>(old_faceCtrlPts, old_rgn_f2);
  auto c011 = get_vector<3>(old_faceCtrlPts, old_rgn_f3);

  LO e0_flip = -1;
  LO e1_flip = -1;
  LO e2_flip = -1;
  LO e3_flip = -1;
  LO e4_flip = -1;
  LO e5_flip = -1;
  {
    auto e0v0 = old_ev2v[old_rgn_e0*2 + 0];
    auto e0v1 = old_ev2v[old_rgn_e0*2 + 1];
    e0_flip = edge_is_flip(e0v0, e0v1, old_rgn_v0, old_rgn_v1);
    auto e1v0 = old_ev2v[old_rgn_e1*2 + 0];
    auto e1v1 = old_ev2v[old_rgn_e1*2 + 1];
    e1_flip = edge_is_flip(e1v0, e1v1, old_rgn_v1, old_rgn_v2);
    auto e2v0 = old_ev2v[old_rgn_e2*2 + 0];
    auto e2v1 = old_ev2v[old_rgn_e2*2 + 1];
    e2_flip = edge_is_flip(e2v0, e2v1, old_rgn_v2, old_rgn_v0);
    auto e3v0 = old_ev2v[old_rgn_e3*2 + 0];
    auto e3v1 = old_ev2v[old_rgn_e3*2 + 1];
    e3_flip = edge_is_flip(e3v0, e3v1, old_rgn_v0, old_rgn_v3);
    auto e4v0 = old_ev2v[old_rgn_e4*2 + 0];
    auto e4v1 = old_ev2v[old_rgn_e4*2 + 1];
    e4_flip = edge_is_flip(e4v0, e4v1, old_rgn_v1, old_rgn_v3);
    auto e5v0 = old_ev2v[old_rgn_e5*2 + 0];
    auto e5v1 = old_ev2v[old_rgn_e5*2 + 1];
    e5_flip = edge_is_flip(e5v0, e5v1, old_rgn_v2, old_rgn_v3);
  }

  auto pts_per_edge = n_edge_pts;
  Real cx100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 0];
  Real cy100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 1];
  Real cz100 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + 2];
  Real cx200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz200 = old_edgeCtrlPts[old_rgn_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e0_flip > 0) {
    swap2(cx100, cx200);
    swap2(cy100, cy200);
    swap2(cz100, cz200);
  }
  Real cx210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 0];
  Real cy210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 1];
  Real cz210 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + 2];
  Real cx120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz120 = old_edgeCtrlPts[old_rgn_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e1_flip > 0) {
    swap2(cx210, cx120);
    swap2(cy210, cy120);
    swap2(cz210, cz120);
  }
  Real cx020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 0];
  Real cy020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 1];
  Real cz020 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + 2];
  Real cx010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz010 = old_edgeCtrlPts[old_rgn_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e2_flip > 0) {
    swap2(cx020, cx010);
    swap2(cy020, cy010);
    swap2(cz020, cz010);
  }
  Real cx001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 0];
  Real cy001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 1];
  Real cz001 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + 2];
  Real cx002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz002 = old_edgeCtrlPts[old_rgn_e3*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e3_flip > 0) {
    swap2(cx002, cx001);
    swap2(cy002, cy001);
    swap2(cz002, cz001);
  }
  Real cx201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 0];
  Real cy201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 1];
  Real cz201 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + 2];
  Real cx102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz102 = old_edgeCtrlPts[old_rgn_e4*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e4_flip > 0) {
    swap2(cx102, cx201);
    swap2(cy102, cy201);
    swap2(cz102, cz201);
  }
  Real cx021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 0];
  Real cy021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 1];
  Real cz021 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + 2];
  Real cx012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz012 = old_edgeCtrlPts[old_rgn_e5*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e5_flip > 0) {
    swap2(cx012, cx021);
    swap2(cy012, cy021);
    swap2(cz012, cz021);
  }
  auto c100 = vector_3(cx100, cy100, cz100);
  auto c200 = vector_3(cx200, cy200, cz200);
  auto c210 = vector_3(cx210, cy210, cz210);
  auto c120 = vector_3(cx120, cy120, cz120);
  auto c020 = vector_3(cx020, cy020, cz020);
  auto c010 = vector_3(cx010, cy010, cz010);
  auto c001 = vector_3(cx001, cy001, cz001);
  auto c002 = vector_3(cx002, cy002, cz002);
  auto c102 = vector_3(cx102, cy102, cz102);
  auto c201 = vector_3(cx201, cy201, cz201);
  auto c012 = vector_3(cx012, cy012, cz012);
  auto c021 = vector_3(cx021, cy021, cz021);

  Vector<3> p11;
  for (LO j = 0; j < dim; ++j) {
    p11[j] = c000[j]*Bijk(order,0,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c300[j]*Bijk(order,3,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c030[j]*Bijk(order,0,3,0,nodePt[0],nodePt[1],nodePt[2]) +
      c003[j]*Bijk(order,0,0,3,nodePt[0],nodePt[1],nodePt[2]) +
      c100[j]*Bijk(order,1,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c200[j]*Bijk(order,2,0,0,nodePt[0],nodePt[1],nodePt[2]) +
      c210[j]*Bijk(order,2,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c120[j]*Bijk(order,1,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c020[j]*Bijk(order,0,2,0,nodePt[0],nodePt[1],nodePt[2]) +
      c010[j]*Bijk(order,0,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c001[j]*Bijk(order,0,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c002[j]*Bijk(order,0,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c201[j]*Bijk(order,2,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c102[j]*Bijk(order,1,0,2,nodePt[0],nodePt[1],nodePt[2]) +
      c012[j]*Bijk(order,0,1,2,nodePt[0],nodePt[1],nodePt[2]) +
      c021[j]*Bijk(order,0,2,1,nodePt[0],nodePt[1],nodePt[2]) +
      c110[j]*Bijk(order,1,1,0,nodePt[0],nodePt[1],nodePt[2]) +
      c101[j]*Bijk(order,1,0,1,nodePt[0],nodePt[1],nodePt[2]) +
      c111[j]*Bijk(order,1,1,1,nodePt[0],nodePt[1],nodePt[2]) +
      c011[j]*Bijk(order,0,1,1,nodePt[0],nodePt[1],nodePt[2]);
  }
  return p11;
}

OMEGA_H_DEVICE Vector<2> face_parametricToParent_2d(
    LO const order, LO const old_face, LOs old_ev2v, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Real nodePts_0, Real nodePts_1, LOs old_fv2v) {
  LO const dim = 2;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e0_flip > 0) {
    swap2(cx10, cx20);
    swap2(cy10, cy20);
  }

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e1_flip > 0) {
    swap2(cx12, cx21);
    swap2(cy12, cy21);
  }

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (e2_flip > 0) {
    swap2(cx02, cx01);
    swap2(cy02, cy01);
  }

  Real cx11 = old_faceCtrlPts[old_face*dim + 0];
  Real cy11 = old_faceCtrlPts[old_face*dim + 1];
  auto c00 = vector_2(cx00, cy00);
  auto c10 = vector_2(cx10, cy10);
  auto c20 = vector_2(cx20, cy20);
  auto c30 = vector_2(cx30, cy30);
  auto c21 = vector_2(cx21, cy21);
  auto c12 = vector_2(cx12, cy12);
  auto c03 = vector_2(cx03, cy03);
  auto c02 = vector_2(cx02, cy02);
  auto c01 = vector_2(cx01, cy01);
  auto c11 = vector_2(cx11, cy11);

  Vector<2> p11_w;
  for (LO k = 0; k < dim; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c30[k]*Bij(order, 3, 0, nodePts_0, nodePts_1) +
      c21[k]*Bij(order, 2, 1, nodePts_0, nodePts_1) +
      c12[k]*Bij(order, 1, 2, nodePts_0, nodePts_1) +
      c03[k]*Bij(order, 0, 3, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

OMEGA_H_DEVICE Vector<3> face_parametricToParent_3d(
    LO const order, LO const old_face, LOs old_ev2v, LOs old_fe2e,
    Reals old_vertCtrlPts, Reals old_edgeCtrlPts, Reals old_faceCtrlPts,
    Real nodePts_0, Real nodePts_1, LOs old_fv2v) {
  LO const dim = 3;
  LO const n_edge_pts = n_internal_ctrlPts(EDGE, order);
  LO const v0_old_face = old_fv2v[old_face*3 + 0];
  LO const v1_old_face = old_fv2v[old_face*3 + 1];
  LO const v2_old_face = old_fv2v[old_face*3 + 2];
  LO const old_face_e0 = old_fe2e[old_face*3 + 0];
  LO const old_face_e1 = old_fe2e[old_face*3 + 1];
  LO const old_face_e2 = old_fe2e[old_face*3 + 2];
  I8 e0_flip = -1;
  I8 e1_flip = -1;
  I8 e2_flip = -1;
  LO v1 = v1_old_face;
  LO v2 = v2_old_face;
  auto e0v0_old_face = old_ev2v[old_face_e0*2 + 0];
  auto e0v1 = old_ev2v[old_face_e0*2 + 1];
  auto e1v0_old_face = old_ev2v[old_face_e1*2 + 0];
  auto e1v1 = old_ev2v[old_face_e1*2 + 1];
  auto e2v0_old_face = old_ev2v[old_face_e2*2 + 0];
  auto e2v1 = old_ev2v[old_face_e2*2 + 1];
  if ((e0v0_old_face == v1) && (e0v1 == v0_old_face)) {
    e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0_old_face == v0_old_face) && (e0v1 == v1));
  }
  if ((e1v0_old_face == v2) && (e1v1 == v1)) {
    e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e1v0_old_face == v1) && (e1v1 == v2));
  }
  if ((e2v0_old_face == v0_old_face) && (e2v1 == v2)) {
    e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e2v0_old_face == v2) && (e2v1 == v0_old_face));
  }

  Real cx00 = old_vertCtrlPts[v0_old_face*dim + 0];
  Real cy00 = old_vertCtrlPts[v0_old_face*dim + 1];
  Real cz00 = old_vertCtrlPts[v0_old_face*dim + 2];
  Real cx30 = old_vertCtrlPts[v1*dim + 0];
  Real cy30 = old_vertCtrlPts[v1*dim + 1];
  Real cz30 = old_vertCtrlPts[v1*dim + 2];
  Real cx03 = old_vertCtrlPts[v2*dim + 0];
  Real cy03 = old_vertCtrlPts[v2*dim + 1];
  Real cz03 = old_vertCtrlPts[v2*dim + 2];

  auto pts_per_edge = n_edge_pts;
  Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
  Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
  Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
  Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e0_flip > 0) {
    auto tempx = cx10;
    auto tempy = cy10;
    auto tempz = cz10;
    cx10 = cx20;
    cy10 = cy20;
    cz10 = cz20;
    cx20 = tempx;
    cy20 = tempy;
    cz20 = tempz;
  }

  Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
  Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
  Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
  Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e1_flip > 0) {
    auto tempx = cx21;
    auto tempy = cy21;
    auto tempz = cz21;
    cx21 = cx12;
    cy21 = cy12;
    cz21 = cz12;
    cx12 = tempx;
    cy12 = tempy;
    cz12 = tempz;
  }

  Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
  Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
  Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
  Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (e2_flip > 0) {
    auto tempx = cx02;
    auto tempy = cy02;
    auto tempz = cz02;
    cx02 = cx01;
    cy02 = cy01;
    cz02 = cz01;
    cx01 = tempx;
    cy01 = tempy;
    cz01 = tempz;
  }

  Real cx11 = old_faceCtrlPts[old_face*dim + 0];
  Real cy11 = old_faceCtrlPts[old_face*dim + 1];
  Real cz11 = old_faceCtrlPts[old_face*dim + 2];
  auto c00 = vector_3(cx00, cy00, cz00);
  auto c10 = vector_3(cx10, cy10, cz10);
  auto c20 = vector_3(cx20, cy20, cz20);
  auto c30 = vector_3(cx30, cy30, cz30);
  auto c21 = vector_3(cx21, cy21, cz21);
  auto c12 = vector_3(cx12, cy12, cz12);
  auto c03 = vector_3(cx03, cy03, cz03);
  auto c02 = vector_3(cx02, cy02, cz02);
  auto c01 = vector_3(cx01, cy01, cz01);
  auto c11 = vector_3(cx11, cy11, cz11);

  Vector<3> p11_w;
  for (LO k = 0; k < 3; ++k) {
    p11_w[k] = c00[k]*Bij(order, 0, 0, nodePts_0, nodePts_1) +
      c10[k]*Bij(order, 1, 0, nodePts_0, nodePts_1) +
      c20[k]*Bij(order, 2, 0, nodePts_0, nodePts_1) +
      c30[k]*Bij(order, 3, 0, nodePts_0, nodePts_1) +
      c21[k]*Bij(order, 2, 1, nodePts_0, nodePts_1) +
      c12[k]*Bij(order, 1, 2, nodePts_0, nodePts_1) +
      c03[k]*Bij(order, 0, 3, nodePts_0, nodePts_1) +
      c02[k]*Bij(order, 0, 2, nodePts_0, nodePts_1) +
      c01[k]*Bij(order, 0, 1, nodePts_0, nodePts_1) +
      c11[k]*Bij(order, 1, 1, nodePts_0, nodePts_1);
  }
  return p11_w;
}

OMEGA_H_DEVICE Vector<2> face_interpToCtrlPt_2d(
    LO const order, LO const newface, LOs new_ev2v, LOs new_fe2e,
    Reals new_vertCtrlPts, Reals new_edgeCtrlPts, Vector<2> p11, LOs new_fv2v) {
  LO const dim=2;
  LO const pts_per_edge = n_internal_ctrlPts(EDGE, order);
  I8 newface_e0_flip = -1;
  I8 newface_e1_flip = -1;
  I8 newface_e2_flip = -1;
  LO newface_v0 = new_fv2v[newface*3 + 0];
  LO newface_v1 = new_fv2v[newface*3 + 1];
  LO newface_v2 = new_fv2v[newface*3 + 2];
  LO newface_e0 = new_fe2e[newface*3 + 0];
  LO newface_e1 = new_fe2e[newface*3 + 1];
  LO newface_e2 = new_fe2e[newface*3 + 2];
  auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
  auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
  auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
  auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
  auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
  auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
  if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
    newface_e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
  }
  if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
    newface_e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
  }
  if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
    newface_e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
  }

  Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
  Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
  Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
  Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
  Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
  Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];

  Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
  Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
  Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e0_flip > 0) {
    swap2(newface_cx10, newface_cx20);
    swap2(newface_cy10, newface_cy20);
  }

  Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
  Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
  Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e1_flip > 0) {
    swap2(newface_cx21, newface_cx12);
    swap2(newface_cy21, newface_cy12);
  }

  Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
  Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
  Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  if (newface_e2_flip > 0) {
    swap2(newface_cx01, newface_cx02);
    swap2(newface_cy01, newface_cy02);
  }
  auto newface_c00 = vector_2(newface_cx00, newface_cy00);
  auto newface_c10 = vector_2(newface_cx10, newface_cy10);
  auto newface_c20 = vector_2(newface_cx20, newface_cy20);
  auto newface_c30 = vector_2(newface_cx30, newface_cy30);
  auto newface_c21 = vector_2(newface_cx21, newface_cy21);
  auto newface_c12 = vector_2(newface_cx12, newface_cy12);
  auto newface_c03 = vector_2(newface_cx03, newface_cy03);
  auto newface_c02 = vector_2(newface_cx02, newface_cy02);
  auto newface_c01 = vector_2(newface_cx01, newface_cy01);

  auto xi_11 = xi_11_cube();
  Vector<2> newface_c11_w;
  for (LO k = 0; k < dim; ++k) {
    newface_c11_w[k] = (p11[k] -
        newface_c00[k]*Bij(order, 0, 0, xi_11[0], xi_11[1]) -
        newface_c10[k]*Bij(order, 1, 0, xi_11[0], xi_11[1]) -
        newface_c20[k]*Bij(order, 2, 0, xi_11[0], xi_11[1]) -
        newface_c30[k]*Bij(order, 3, 0, xi_11[0], xi_11[1]) -
        newface_c21[k]*Bij(order, 2, 1, xi_11[0], xi_11[1]) -
        newface_c12[k]*Bij(order, 1, 2, xi_11[0], xi_11[1]) -
        newface_c03[k]*Bij(order, 0, 3, xi_11[0], xi_11[1]) -
        newface_c02[k]*Bij(order, 0, 2, xi_11[0], xi_11[1]) -
        newface_c01[k]*Bij(order, 0, 1, xi_11[0], xi_11[1]))/
      Bij(order, 1, 1, xi_11[0], xi_11[1]);
  }

  return newface_c11_w;
}

OMEGA_H_DEVICE Vector<3> face_interpToCtrlPt_3d(
    LO const order, LO const newface, LOs new_ev2v, LOs new_fe2e,
    Reals new_vertCtrlPts, Reals new_edgeCtrlPts, Vector<3> p11, LOs new_fv2v) {
  LO const dim=3;
  LO const pts_per_edge = n_internal_ctrlPts(EDGE, order);
  I8 newface_e0_flip = -1;
  I8 newface_e1_flip = -1;
  I8 newface_e2_flip = -1;
  LO newface_v0 = new_fv2v[newface*3 + 0];
  LO newface_v1 = new_fv2v[newface*3 + 1];
  LO newface_v2 = new_fv2v[newface*3 + 2];
  LO newface_e0 = new_fe2e[newface*3 + 0];
  LO newface_e1 = new_fe2e[newface*3 + 1];
  LO newface_e2 = new_fe2e[newface*3 + 2];
  auto newface_e0v0 = new_ev2v[newface_e0*2 + 0];
  auto newface_e0v1 = new_ev2v[newface_e0*2 + 1];
  auto newface_e1v0 = new_ev2v[newface_e1*2 + 0];
  auto newface_e1v1 = new_ev2v[newface_e1*2 + 1];
  auto newface_e2v0 = new_ev2v[newface_e2*2 + 0];
  auto newface_e2v1 = new_ev2v[newface_e2*2 + 1];
  if ((newface_e0v0 == newface_v1) && (newface_e0v1 == newface_v0)) {
    newface_e0_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e0v0 == newface_v0) && (newface_e0v1 == newface_v1));
  }
  if ((newface_e1v0 == newface_v2) && (newface_e1v1 == newface_v1)) {
    newface_e1_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e1v0 == newface_v1) && (newface_e1v1 == newface_v2));
  }
  if ((newface_e2v0 == newface_v0) && (newface_e2v1 == newface_v2)) {
    newface_e2_flip = 1;
  }
  else {
    OMEGA_H_CHECK((newface_e2v0 == newface_v2) && (newface_e2v1 == newface_v0));
  }

  Real newface_cx00 = new_vertCtrlPts[newface_v0*dim + 0];
  Real newface_cy00 = new_vertCtrlPts[newface_v0*dim + 1];
  Real newface_cz00 = new_vertCtrlPts[newface_v0*dim + 2];
  Real newface_cx30 = new_vertCtrlPts[newface_v1*dim + 0];
  Real newface_cy30 = new_vertCtrlPts[newface_v1*dim + 1];
  Real newface_cz30 = new_vertCtrlPts[newface_v1*dim + 2];
  Real newface_cx03 = new_vertCtrlPts[newface_v2*dim + 0];
  Real newface_cy03 = new_vertCtrlPts[newface_v2*dim + 1];
  Real newface_cz03 = new_vertCtrlPts[newface_v2*dim + 2];

  Real newface_cx10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 0];
  Real newface_cy10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 1];
  Real newface_cz10 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + 2];
  Real newface_cx20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz20 = new_edgeCtrlPts[newface_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e0_flip > 0) {
    swap2(newface_cx10, newface_cx20);
    swap2(newface_cy10, newface_cy20);
    swap2(newface_cz10, newface_cz20);
  }

  Real newface_cx21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 0];
  Real newface_cy21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 1];
  Real newface_cz21 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + 2];
  Real newface_cx12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz12 = new_edgeCtrlPts[newface_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e1_flip > 0) {
    swap2(newface_cx21, newface_cx12);
    swap2(newface_cy21, newface_cy12);
    swap2(newface_cz21, newface_cz12);
  }

  Real newface_cx02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 0];
  Real newface_cy02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 1];
  Real newface_cz02 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + 2];
  Real newface_cx01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
  Real newface_cy01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
  Real newface_cz01 = new_edgeCtrlPts[newface_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
  if (newface_e2_flip > 0) {
    swap2(newface_cx01, newface_cx02);
    swap2(newface_cy01, newface_cy02);
    swap2(newface_cz01, newface_cz02);
  }
  auto newface_c00 = vector_3(newface_cx00, newface_cy00, newface_cz00);
  auto newface_c10 = vector_3(newface_cx10, newface_cy10, newface_cz10);
  auto newface_c20 = vector_3(newface_cx20, newface_cy20, newface_cz20);
  auto newface_c30 = vector_3(newface_cx30, newface_cy30, newface_cz30);
  auto newface_c21 = vector_3(newface_cx21, newface_cy21, newface_cz21);
  auto newface_c12 = vector_3(newface_cx12, newface_cy12, newface_cz12);
  auto newface_c03 = vector_3(newface_cx03, newface_cy03, newface_cz03);
  auto newface_c02 = vector_3(newface_cx02, newface_cy02, newface_cz02);
  auto newface_c01 = vector_3(newface_cx01, newface_cy01, newface_cz01);

  auto xi_11 = xi_11_cube();
  Vector<3> newface_c11_w;
  for (LO k = 0; k < 3; ++k) {
    newface_c11_w[k] = (p11[k] -
        newface_c00[k]*Bij(order, 0, 0, xi_11[0], xi_11[1]) -
        newface_c10[k]*Bij(order, 1, 0, xi_11[0], xi_11[1]) -
        newface_c20[k]*Bij(order, 2, 0, xi_11[0], xi_11[1]) -
        newface_c30[k]*Bij(order, 3, 0, xi_11[0], xi_11[1]) -
        newface_c21[k]*Bij(order, 2, 1, xi_11[0], xi_11[1]) -
        newface_c12[k]*Bij(order, 1, 2, xi_11[0], xi_11[1]) -
        newface_c03[k]*Bij(order, 0, 3, xi_11[0], xi_11[1]) -
        newface_c02[k]*Bij(order, 0, 2, xi_11[0], xi_11[1]) -
        newface_c01[k]*Bij(order, 0, 1, xi_11[0], xi_11[1]))/
      Bij(order, 1, 1, xi_11[0], xi_11[1]);

  }

  return newface_c11_w;
}

LOs create_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges);

void create_curved_faces_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

LOs create_curved_verts_and_edges_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts,
                                     LOs keys2edges);

void create_curved_faces_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

LOs coarsen_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh,
                                      LOs old_ents2new_ents, LOs prods2new_ents,
                                      LOs keys2prods, LOs old_verts2new_verts, LOs old_edges2new_edges);

template <Int ent_dim>
LOs checkValidity(Mesh *new_mesh, LOs new_tris);

#define OMEGA_H_EXPL_INST_DECL(ent_dim)                                            \
  extern template LOs checkValidity(Mesh *new_mesh, LOs new_tris);
OMEGA_H_EXPL_INST_DECL(1)
OMEGA_H_EXPL_INST_DECL(2)
OMEGA_H_EXPL_INST_DECL(3)
#undef OMEGA_H_EXPL_INST_DECL

} // namespace Omega_h

#endif
