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

OMEGA_H_DEVICE cubic_region_xi_values
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
  else if (old_key_edge == old_e1) {
    p11[0] = 1.0/3.0;
    p11[1] = 1.0/6.0;
    p11[2] = 1.0/6.0;
  }
  else {
    Omega_h_fail("incorrect old rgn\n");
  }

  return Reals(p11);
}

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
