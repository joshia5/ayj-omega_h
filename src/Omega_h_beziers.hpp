#ifndef OMEGA_H_BEZIERS_HPP
#define OMEGA_H_BEZIERS_HPP

#include "Omega_h_array.hpp"

namespace Omega_h {

Real B0_quad(Real u);
Real B1_quad(Real u);
Real B2_quad(Real u);

Real B0_cube(Real u);
Real B1_cube(Real u);
Real B2_cube(Real u);
Real B3_cube(Real u);

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

//TODO this fn will need to be called from both host and device
constexpr OMEGA_H_DEVICE Real 
Bijk(LO const P, LO const i, LO const j, LO const k, Real const u,
     Real const v, Real const w) noexcept {
  LO const l = P - i - j - k;
  OMEGA_H_CHECK(l >= 0);
  Real const t = 1.0 - u - v - w;
  OMEGA_H_CHECK((t >= 0.0) && (t <= 1.0));
  return factorial(1.0*P)*
    std::pow(u,i)*std::pow(v,j)*std::pow(w,k)*std::pow(t,l)/(
    factorial(1.0*i)*factorial(1.0*j)*factorial(1.0*k)*factorial(1.0*l));
}

constexpr OMEGA_H_DEVICE Real 
Bij(LO const P, LO const i, LO const j, Real const u, Real const v) noexcept {
  LO const k = P - i - j;
  OMEGA_H_CHECK((k >= 0) && (k <= P));
  Real const w = 1.0 - u - v;
  OMEGA_H_CHECK((w >= 0.0) && (w <= 1.0));
  return factorial(1.0*P)*
    std::pow(u,i)*std::pow(v,j)*std::pow(w,k)/(
    factorial(1.0*i)*factorial(1.0*j)*factorial(1.0*k));
}

constexpr OMEGA_H_DEVICE Real 
Bi(LO const P, LO const i, Real const u) noexcept {
  LO const j = P - i;
  OMEGA_H_CHECK((j >= 0) && (j <= P));
  Real const v = 1.0 - u;
  OMEGA_H_CHECK((v >= 0.0) && (v <= 1.0));
  return factorial(1.0*P)*
    std::pow(u,i)*std::pow(v,j)/(factorial(1.0*i)*factorial(1.0*j));
}

Real xi_1_quad();

Real xi_1_cube();
Real xi_2_cube();

Real xi_1_quart();
Real xi_2_quart();
Real xi_3_quart();

Real xi_1_quint();
Real xi_2_quint();
Real xi_3_quint();
Real xi_4_quint();

Reals xi_11_cube();

void elevate_curve_order_2to3(Mesh* mesh);
void elevate_curve_order_3to4(Mesh* mesh);
void elevate_curve_order_4to5(Mesh* mesh);
void elevate_curve_order_5to6(Mesh* mesh);

void calc_quad_ctrlPts_from_interpPts(Mesh *mesh);

OMEGA_H_DEVICE Reals cubic_noKeyEdge_xi_values(LO old_vert, LO v0, LO v1, LO v2,
                                        LO old_edge, LO e0, LO e1, LO e2) {

  Write<Real> p1_p2(4);
  if (old_vert == v0) {
    OMEGA_H_CHECK(old_edge == e1);
    p1_p2[0] = 0.13740215;
    p1_p2[1] = 0.13740215;
    p1_p2[2] = 0.362597849;
    p1_p2[3] = 0.362597849;
  }
  else if (old_vert == v1) {
    OMEGA_H_CHECK(old_edge == e2);
    p1_p2[0] = 0.7251957;
    p1_p2[1] = 0.13740215;
    p1_p2[2] = 0.2748043;
    p1_p2[3] = 0.362597849;
  }
  else if (old_vert == v2) {
    OMEGA_H_CHECK(old_edge == e0);
    p1_p2[0] = 0.13740215;
    p1_p2[1] = 0.7251957;
    p1_p2[2] = 0.362597849;
    p1_p2[3] = 0.2748043;
  }
  else {
    Omega_h_fail("incorrect old face\n");
  }

 return Reals(p1_p2);
}

OMEGA_H_DEVICE Reals cubic_face_xi_values
  (LO const old_noKey_vert, LO const old_v0, LO const old_v1, LO const old_v2,
   LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const new_v0_f0, LO const new_v1_f0, LO const new_v2_f0,
   LO const new_v0_f1, LO const new_v1_f1, LO const new_v2_f1,
   LO const oldf_v0_new, LO const oldf_v1_new, LO const oldf_v2_new) {

  Write<Real> p1_p2(4);
  I8 should_swap = -1;
  if (old_noKey_vert == old_v0) {
    OMEGA_H_CHECK(old_key_edge == old_e1);
    p1_p2[0] = 1.0/2.0;
    p1_p2[1] = 1.0/6.0;
    p1_p2[2] = 1.0/6.0;
    p1_p2[3] = 1.0/2.0;

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
    p1_p2[0] = 1.0/3.0;
    p1_p2[1] = 1.0/2.0;
    p1_p2[2] = 1.0/3.0;
    p1_p2[3] = 1.0/6.0;

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
    p1_p2[0] = 1.0/6.0;
    p1_p2[1] = 1.0/3.0;
    p1_p2[2] = 1.0/2.0;
    p1_p2[3] = 1.0/3.0;

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
    Omega_h_fail("incorrect old face\n");
  }
  if (should_swap == 1) {
    swap2(p1_p2[0], p1_p2[2]);
    swap2(p1_p2[1], p1_p2[3]);
  }

  return Reals(p1_p2);
}

OMEGA_H_DEVICE Reals cubic_region_xi_values
  (LO const old_key_edge, LO const old_e0, LO const old_e1, LO const old_e2,
   LO const old_e3, LO const old_e4, LO const old_e5) {
  Write<Real> p11(3);
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
    Omega_h_fail("incorrect old rgn\n");
  }

  return Reals(p11);
}

OMEGA_H_DEVICE LO edge_is_flip(LO const e0v0, LO const e0v1, LO const v0, LO const v1) {
  LO is_flip = -1;
  if ((e0v0 == v1) && (e0v1 == v0)) {
    is_flip = 1;
  }
  else {
    OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
  }
  return is_flip;
}

OMEGA_H_DEVICE Reals face_interpToCtrlPt_3d(LO const order, LO const newface, LOs new_ev2v, LOs new_fe2e,
    Reals new_vertCtrlPts, Reals new_edgeCtrlPts, Reals p11, LOs new_fv2v) {
  LO const dim=3;
  LO const pts_per_edge=2;
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
  Write<Real> newface_c11_w(3);
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

    //TODO mult pts per face
  }

  return Reals(newface_c11_w);
}

//Reals OMEGA_H_DEVICE coordsFromXi(Int ent_dim, 

LOs create_curved_verts_and_edges_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts);

void create_curved_faces_2d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

LOs create_curved_verts_and_edges_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new,
                                     LOs prods2new, LOs keys2prods,
                                     LOs keys2midverts, LOs old_verts2new_verts);

void create_curved_faces_3d(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                            LOs keys2prods, LOs keys2edges, LOs keys2old_faces,
                            LOs old_verts2new_verts);

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
