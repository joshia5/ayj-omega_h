#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_vector.hpp"

namespace Omega_h {

Real B0_quad(Real u) {
  return (1.0-u)*(1.0-u);
}

Real B1_quad(Real u) {
  return 2.0*u*(1.0-u);
}

Real B2_quad(Real u) {
  return u*u;
}

Real B0_cube(Real u) {
  return (1.0-u)*(1.0-u)*(1.0-u);
}

Real B1_cube(Real u) {
  return 3.0*u*(1.0-u)*(1.0-u);
}

Real B2_cube(Real u) {
  return 3.0*(1.0-u)*u*u;
}

Real B3_cube(Real u) {
  return u*u*u;
}

Real B0_quart(Real u) {
  return (1.0-u)*(1.0-u)*(1.0-u)*(1.0-u);
}

Real B1_quart(Real u) {
  return 4.0*u*(1.0-u)*(1.0-u)*(1.0-u);
}

Real B2_quart(Real u) {
  return 6.0*u*u*(1.0-u)*(1.0-u);
}

Real B3_quart(Real u) {
  return 4.0*u*u*u*(1.0-u);
}

Real B4_quart(Real u) {
  return u*u*u*u;
}

Real B00_quad(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v);
}

Real B10_quad(Real u, Real v) {
  return 2.0*u*(1.0-u-v);
}

Real B20_quad(Real u, Real v) {
  return u*u + 0.0*v;
}

Real B11_quad(Real u, Real v) {
  return 2.0*u*v;
}

Real B02_quad(Real u, Real v) {
  return v*v + 0.0*u;
}

Real B01_quad(Real u, Real v) {
  return 2.0*v*(1.0-u-v);
}

Real B00_cube(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B10_cube(Real u, Real v) {
  return 3.0*u*(1.0-u-v)*(1.0-u-v);
}

Real B20_cube(Real u, Real v) {
  return 3.0*u*u*(1.0-u-v);
}

Real B30_cube(Real u, Real v) {
  return u*u*u + 0.0*v;
}

Real B21_cube(Real u, Real v) {
  return 3.0*u*u*v;
}

Real B12_cube(Real u, Real v) {
  return 3.0*u*v*v;
}

Real B03_cube(Real u, Real v) {
  return v*v*v + 0.0*u;
}

Real B02_cube(Real u, Real v) {
  return 3.0*v*v*(1.0-u-v);
}

Real B01_cube(Real u, Real v) {
  return 3.0*v*(1.0-u-v)*(1.0-u-v);
}

Real B11_cube(Real u, Real v) {
  return 6.0*u*v*(1.0-u-v);
}

Real B00_quart(Real u, Real v) {
  return (1.0-u-v)*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B10_quart(Real u, Real v) {
  return 4.0*u*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B20_quart(Real u, Real v) {
  return 6.0*u*u*(1.0-u-v)*(1.0-u-v);
}

Real B30_quart(Real u, Real v) {
  return 4.0*u*u*u*(1.0-u-v);
}

Real B40_quart(Real u, Real v) {
  return u*u*u*u + 0.0*v;
}

Real B31_quart(Real u, Real v) {
  return 4.0*u*u*u*v;
}

Real B22_quart(Real u, Real v) {
  return 6.0*u*u*v*v;
}

Real B13_quart(Real u, Real v) {
  return 4.0*u*v*v*v;
}

Real B04_quart(Real u, Real v) {
  return v*v*v*v + 0.0*u;
}

Real B03_quart(Real u, Real v) {
  return 4.0*v*v*v*(1.0-u-v);
}

Real B02_quart(Real u, Real v) {
  return 6.0*v*v*(1.0-u-v)*(1.0-u-v);
}

Real B01_quart(Real u, Real v) {
  return 4.0*v*(1.0-u-v)*(1.0-u-v)*(1.0-u-v);
}

Real B11_quart(Real u, Real v) {
  return 12.0*u*v*(1.0-u-v)*(1.0-u-v);
}

Real B21_quart(Real u, Real v) {
  return 12.0*u*u*v*(1.0-u-v);
}

Real B12_quart(Real u, Real v) {
  return 12.0*u*v*v*(1.0-u-v);
}

Real xi_1_quad() {
 return 0.5;
}

// Babushka chen points
Real xi_1_cube() {
 return 0.2748043;
}
Real xi_2_cube() {
  return 0.7251957;
}

Real xi_1_quart() {
 return 0.1693976;
}
Real xi_2_quart() {
 return 0.5;
}
Real xi_3_quart() {
 return 0.8306024;
}

Real xi_1_quint() {
 return 0.1133573;
}
Real xi_2_quint() {
 return 0.3568239;
}
Real xi_3_quint() {
 return 0.6431761;
}
Real xi_4_quint() {
  return 0.8866427;
}

OMEGA_H_DEVICE Reals xi_values_cubic(LO old_vert, LO v0, LO v1, LO v2) {

  Write<Real> p1_p2_p3(6);
  if (old_vert == v0) {
    p1_p2_p3[0] = 0.13740215;
    p1_p2_p3[1] = 0.13740215;
    p1_p2_p3[2] = 0.362597849;
    p1_p2_p3[3] = 0.362597849;

    p1_p2_p3[4] = 0.5;
    p1_p2_p3[5]= 0.5;
  }
  else if (old_vert == v1) {
    p1_p2_p3[2] = 0.2748043;
    p1_p2_p3[3] = 0.362597849;
    p1_p2_p3[0] = 0.7251957;
    p1_p2_p3[1] = 0.13740215;

    p1_p2_p3[4] = 0.0;
    p1_p2_p3[5]= 0.5;
  }
  else  {
    OMEGA_H_CHECK(old_vert == v2);
    p1_p2_p3[0] = 0.7251957;
    p1_p2_p3[1] = 0.13740215;
    p1_p2_p3[2] = 0.2748043;
    p1_p2_p3[3] = 0.362597849;

    p1_p2_p3[4] = 0.5;
    p1_p2_p3[5]= 0.0;
  }

  return Reals(p1_p2_p3);
}

void calc_quad_ctrlPts_from_interpPts(Mesh *mesh) {
  auto coords = mesh->coords();
  auto interpPts = mesh->get_ctrlPts(1);
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto dim = mesh->dim();
  auto nedge = mesh->nedges();
  Real xi_1 = 0.5;
  Write<Real> new_pts(nedge*dim, 0.0);

  auto f = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];

    //TODO make this a templated fn like cal_pts_tmpl
    if (dim == 2) {
      Vector<2> c1;
      auto c0 = get_vector<2>(coords, v0);
      auto p1 = get_vector<2>(interpPts, i);
      auto c2 = get_vector<2>(coords, v1);
      for (Int j = 0; j < dim; ++j) {
        c1[j] = (p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c2[j])/B1_quad(xi_1);
      }
      set_vector(new_pts, i, c1);
    }
    else {
      OMEGA_H_CHECK(dim == 3);
      Vector<3> c1;
      auto c0 = get_vector<3>(coords, v0);
      auto p1 = get_vector<3>(interpPts, i);
      auto c2 = get_vector<3>(coords, v1);
      for (Int j = 0; j < dim; ++j) {
        c1[j] = (p1[j] - B0_quad(xi_1)*c0[j] - B2_quad(xi_1)*c2[j])/B1_quad(xi_1);
      }
      set_vector(new_pts, i, c1);
    }
  };
  parallel_for(nedge, std::move(f));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

void elevate_curve_order_2to3(Mesh* mesh) {
  I8 new_order = 3;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  auto calc_edge_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/3.0)*coords[v0*dim + d] +
              (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (1.0/3.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + dim + d] = c2[d];
    }
  };
  parallel_for(nedge, calc_edge_pts);
  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));

  auto nface = mesh->nfaces();
  n_new_pts = mesh->n_internal_ctrlPts(2);
  Write<Real> face_pts(nface*n_new_pts*dim, 0.0);
  Write<Real> c11(dim, 0.0);
  auto calc_face_pts = OMEGA_H_LAMBDA (LO i) {
    auto e0 = fe2e[i*3];
    auto e1 = fe2e[i*3 + 1];
    auto e2 = fe2e[i*3 + 2];
    for (LO d = 0; d < dim; ++d) {
      c11[d] = (1.0/3.0)*old_ctrl_pts[e0*old_n_ctrl_pts*dim + d] +
               (1.0/3.0)*old_ctrl_pts[e1*old_n_ctrl_pts*dim + d] +
               (1.0/3.0)*old_ctrl_pts[e2*old_n_ctrl_pts*dim + d];
      face_pts[i*n_new_pts*dim + d] = c11[d];
    }
  };
  parallel_for(nface, calc_face_pts);
  mesh->set_tag_for_ctrlPts(2, Reals(face_pts));

  return;
}

void elevate_curve_order_3to4(Mesh* mesh) {
  I8 new_order = 4;
  auto old_edge_ctrl_pts = mesh->get_ctrlPts(1);
  auto pts_per_edge = mesh->n_internal_ctrlPts(1);
  auto old_face_ctrlPts = mesh->get_ctrlPts(2);
  auto old_face_n_ctrl_pts = mesh->n_internal_ctrlPts(2);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto fv2v = mesh->ask_down(2, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/4.0)*coords[v0*dim + d] +
              (3.0/4.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + d];
      c2[d] = (1.0/2.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + d] +
              (1.0/2.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + 1*dim + d];
      c3[d] = (3.0/4.0)*old_edge_ctrl_pts[i*pts_per_edge*dim + 1*dim + d] +
              (1.0/4.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));

  auto nface = mesh->nfaces();
  auto n_new_face_pts = mesh->n_internal_ctrlPts(2);
  Write<Real> face_pts(nface*n_new_face_pts*dim, 0.0);

  auto calc_face_pts = OMEGA_H_LAMBDA (LO i) {
    auto e0 = fe2e[i*3];
    auto e1 = fe2e[i*3 + 1];
    auto e2 = fe2e[i*3 + 2];
    auto v0 = fv2v[i*3];
    auto v1 = fv2v[i*3 + 1];
    auto v2 = fv2v[i*3 + 2];

    I8 e0_flip = -1;
    I8 e1_flip = -1;
    I8 e2_flip = -1;

    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    if ((e0v0 == v1) && (e0v1 == v0)) {
      e0_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      e1_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      e2_flip = 1;
    }
    else {
      OMEGA_H_CHECK((e2v0 == v2) && (e2v1 == v0));
    }

    Real cx10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 0];
    Real cy10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 1];
    Real cz10 = old_edge_ctrl_pts[e0*pts_per_edge*dim + 2];
    Real cx20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 0];
    Real cy20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 1];
    Real cz20 = old_edge_ctrl_pts[e0*pts_per_edge*dim + dim + 2];
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

    Real cx21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 0];
    Real cy21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 1];
    Real cz21 = old_edge_ctrl_pts[e1*pts_per_edge*dim + 2];
    Real cx12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 0];
    Real cy12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 1];
    Real cz12 = old_edge_ctrl_pts[e1*pts_per_edge*dim + dim + 2];
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

    Real cx02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 0];
    Real cy02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 1];
    Real cz02 = old_edge_ctrl_pts[e2*pts_per_edge*dim + 2];
    Real cx01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 0];
    Real cy01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 1];
    Real cz01 = old_edge_ctrl_pts[e2*pts_per_edge*dim + dim + 2];
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

    auto oldc11 = Reals({old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 0],
                         old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 1],
                         old_face_ctrlPts[i*old_face_n_ctrl_pts*dim + 2]});

    Write<Real> c11(dim, 0.0);
    Write<Real> c21(dim, 0.0);
    Write<Real> c12(dim, 0.0);

    c11[0] = (1.0/4.0)*(cx10 + cx01 + 2.0*oldc11[0]);
    c11[1] = (1.0/4.0)*(cy10 + cy01 + 2.0*oldc11[1]);
    c11[2] = (1.0/4.0)*(cz10 + cz01 + 2.0*oldc11[2]);
    c21[0] = (1.0/4.0)*(cx12 + cx02 + 2.0*oldc11[0]);
    c21[1] = (1.0/4.0)*(cy12 + cy02 + 2.0*oldc11[1]);
    c21[2] = (1.0/4.0)*(cz12 + cz02 + 2.0*oldc11[2]);
    c12[0] = (1.0/4.0)*(cx20 + cx21 + 2.0*oldc11[0]);
    c12[1] = (1.0/4.0)*(cy20 + cy21 + 2.0*oldc11[1]);
    c12[2] = (1.0/4.0)*(cz20 + cz21 + 2.0*oldc11[2]);

    for (LO d = 0; d < dim; ++d) {
      face_pts[i*n_new_face_pts*dim + d] = c11[d];
      face_pts[i*n_new_face_pts*dim + dim + d] = c21[d];
      face_pts[i*n_new_face_pts*dim + dim + dim + d] = c12[d];
    }
  };
  parallel_for(nface, std::move(calc_face_pts));
  mesh->set_tag_for_ctrlPts(2, Reals(face_pts));

  //formulae to calc pts inside regions?

  return;
}

void elevate_curve_order_4to5(Mesh* mesh) {
  I8 new_order = 5;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  Write<Real> c4(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/5.0)*coords[v0*dim + d] +
              (4.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (2.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (3.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d];
      c3[d] = (3.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d] +
              (2.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d];
      c4[d] = (4.0/5.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d] +
              (1.0/5.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
      new_pts[i*n_new_pts*dim + 3*dim + d] = c4[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

void elevate_curve_order_5to6(Mesh* mesh) {
  I8 new_order = 6;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->change_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  Write<Real> c3(dim, 0.0);
  Write<Real> c4(dim, 0.0);
  Write<Real> c5(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
    auto v0 = ev2v[i*2];
    auto v1 = ev2v[i*2 + 1];
    for (LO d = 0; d < dim; ++d) {
      c1[d] = (1.0/6.0)*coords[v0*dim + d] +
              (5.0/6.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (1.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d];
      c3[d] = (1.0/2.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d] +
              (1.0/2.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d];
      c4[d] = (2.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 2*dim + d] +
              (1.0/3.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 3*dim + d];
      c5[d] = (5.0/6.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 3*dim + d] +
              (1.0/6.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
      new_pts[i*n_new_pts*dim + 3*dim + d] = c4[d];
      new_pts[i*n_new_pts*dim + 4*dim + d] = c5[d];
    }
  };
  parallel_for(nedge, std::move(calc_pts));

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

void create_curved_edges(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                         LOs keys2prods, LOs keys2midverts,
                         LOs old_verts2new_verts) {
  printf("in create curved edges fn\n");
  OMEGA_H_TIME_FUNCTION;
  auto nold_edge = old2new.size();
  auto old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto old_coords = mesh->coords();
  auto old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto order = mesh->get_max_order();
  auto dim = mesh->dim();
  auto n_edge_pts = mesh->n_internal_ctrlPts(1);

  auto new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto new_coords = new_mesh->coords();
  auto nnew_edge = new_mesh->nedges();
  auto nnew_verts = new_mesh->nverts();

  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim);
  Write<Real> vert_ctrlPts(nnew_verts*1*dim);
  auto copy_coords = OMEGA_H_LAMBDA(LO i) {
    vert_ctrlPts[i] = new_coords[i];
  };
  parallel_for(new_coords.size(), copy_coords, "init vtx ctrlPts");
  //TODO faces cant be added until the bottom up algo reaches faces
  //auto n_face_pts = mesh->n_internal_ctrlPts(2);

  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  Write<LO> count_key(1, 0);
  OMEGA_H_CHECK(order == 3);
  //Write<Real> edge_curveXi(nnew_edge*2*1);

  auto create_crv_edges = OMEGA_H_LAMBDA (LO old_edge) {
    LO v0_old = old_ev2v[old_edge*2 + 0];
    LO v1_old = old_ev2v[old_edge*2 + 1];

    if (old2new[old_edge] != -1) {
      // map as is from old id to new id
      LO new_edge = old2new[old_edge];
      printf("old edge %d same as new edge %d\n", old_edge, new_edge);
      for (I8 d = 0; d < dim; ++d) {
        edge_ctrlPts[new_edge*n_edge_pts*dim + d] = 
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + d];
        edge_ctrlPts[new_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d] =
          old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + d];
      }
    }
    else {
      LO key_id = count_key[0];
      LO mid_vert = keys2midverts[key_id];
      LO start = keys2prods[key_id];
      LO end = keys2prods[key_id + 1] - 1;
      OMEGA_H_CHECK((end-start) == 2);
      LO new_e0 = prods2new[start];
      LO new_e1 = prods2new[start+1];
      LO new_e2 = prods2new[start+2];

      LO v1_new_e0 = new_ev2v[new_e0*2 + 1];
      LO v0_new_e1 = new_ev2v[new_e1*2 + 0];
      OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));
      //LO v0_new_e0 = new_ev2v[new_e0*2 + 0];
      //LO v1_new_e1 = new_ev2v[new_e1*2 + 1];
      //OMEGA_H_CHECK((v0_new_e0 == v0_old) && (v1_new_e1 == v1_old));
      //topologically yes but diff mesh objects

      LO v0_new_e2 = new_ev2v[new_e2*2 + 0];
      LO v1_new_e2 = new_ev2v[new_e2*2 + 1];
      OMEGA_H_CHECK(v1_new_e2 == mid_vert);
      // swap unneeded, order checked

      /*
      // check flip
      auto flip_new_e0 = -1;
      auto flip_new_e1 = -1;
      */

      //ctrl pts for e0
      {
        //for 2d mesh for now, order=3
        //query old ctrl pts
        Real cx0 = old_coords[v0_old*dim + 0];
        Real cy0 = old_coords[v0_old*dim + 1];
        Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 0];
        Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 1];
        Real cx3 = old_coords[v1_old*dim + 0];
        Real cy3 = old_coords[v1_old*dim + 1];

        Real new_xi_3 = 0.5;
        Real new_cx3 = cx0*B0_cube(new_xi_3) + cx1*B1_cube(new_xi_3) +
                       cx2*B2_cube(new_xi_3) + cx3*B3_cube(new_xi_3);
        Real new_cy3 = cy0*B0_cube(new_xi_3) + cy1*B1_cube(new_xi_3) +
                       cy2*B2_cube(new_xi_3) + cy3*B3_cube(new_xi_3);
        printf(" new mid interp pt %f, %f \n", new_cx3, new_cy3);
        vert_ctrlPts[mid_vert*1*dim + 0] = new_cx3;
        vert_ctrlPts[mid_vert*1*dim + 1] = new_cy3;
        //vert_ctrlPts[mid_vert*1*dim + 2] = new_cz3;

        Real old_xi_2 = xi_2_cube();
        Real new_xi_2 = old_xi_2/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);
        printf("edge 0 p2 %f, %f , xi2 %f\n", new_px2, new_py2, new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = old_xi_1/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);
        printf("edge 0 p1 %f, %f \n", new_px1, new_py1);

        Matrix<2,1> c0({cx0, cy0});
        Matrix<2,1> c3({new_cx3, new_cy3});

        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        /*
        // Doesnt change output, are matrics symmetric?
        // Yes , as we have a cubic case. trying to think will they always be
        // for even higher order??
        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B1_cube(old_xi_2), B2_cube(old_xi_1),
                            B2_cube(old_xi_2)});
        Matrix<2,2> M2({B0_cube(old_xi_1), B0_cube(old_xi_2), B3_cube(old_xi_1),
                        B3_cube(old_xi_2)});
                        */

        auto M1 = invert(M1_inv);

        Matrix<2,1> cx({cx0, new_cx3});
        Matrix<2,1> cy({cy0, new_cy3});
        
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        //auto Cx = M1*fx - M1*M2*c0;
        //auto Cy = M1*fy - M1*M2*c3;
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        printf("\n edge 0 c1 %f, %f \n", Cx(0,0), Cy(0,0));
        printf(" edge 0 c2 %f, %f\n", Cx(1,0), Cy(1,0));

        //edge_curveXi[new_e0*2 + 0] = 0.0;
        //edge_curveXi[new_e0*2 + 1] = 1.0;
      }
      //ctrl pts for e1
      {
        //for 2d mesh for now, order=3
        //query old ctrl pts
        Real cx0 = old_coords[v0_old*dim + 0];
        Real cy0 = old_coords[v0_old*dim + 1];
        Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0];
        Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1];
        Real cx3 = old_coords[v1_old*dim + 0];
        Real cy3 = old_coords[v1_old*dim + 1];

        Real old_xi_2 = xi_2_cube();
        Real new_xi_2 = 0.5 + old_xi_2/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = 0.5 + old_xi_1/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);

        Matrix<2,1> cx({vert_ctrlPts[mid_vert*1*dim + 0], cx3});
        Matrix<2,1> cy({vert_ctrlPts[mid_vert*1*dim + 1], cy3});
                        
        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        auto M1 = invert(M1_inv);
        auto Cx = M1*fx - M1*M2*cx;
        auto Cy = M1*fy - M1*M2*cy;
        edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        //edge_curveXi[new_e1*2 + 0] = 0.0;
        //edge_curveXi[new_e1*2 + 1] = 1.0;
      }
      //ctrl pts for e2
      {
        auto ab2b = new_verts2old_verts.ab2b;
        auto a2ab = new_verts2old_verts.a2ab;
        LO old_vert_noKey = ab2b[a2ab[v0_new_e2]];

        LO old_face = -1;
        for (LO index = old_e2ef[old_edge]; index < old_e2ef[old_edge + 1];
             ++index) {
          LO adj_face = old_ef2f[index];
          for (LO vert = 0; vert < 3; ++vert) {
            LO vert_old_face = old_fv2v[adj_face*3 + vert];
            if (vert_old_face == old_vert_noKey) {
              old_face = adj_face;
              break;
            }
          }
          if (old_face > 0) {
            break;
          }
        }

        // for that old-face, query points at calculated xi values
        printf("For old edge %d, found old face %d\n", old_edge, old_face);
        {
          auto v0_old_face = old_fv2v[old_face*3];
          auto v1 = old_fv2v[old_face*3 + 1];
          auto v2 = old_fv2v[old_face*3 + 2];

          //TODO check alignment of new noKey edge
          //from the example, for now p1 is nearer to noKey vert
          auto nodePts = xi_values_cubic(old_vert_noKey, v0_old_face, v1, v2);

          auto old_face_e0 = old_fe2e[old_face*3];
          auto old_face_e1 = old_fe2e[old_face*3 + 1];
          auto old_face_e2 = old_fe2e[old_face*3 + 2];

          I8 e0_flip = -1;
          I8 e1_flip = -1;
          I8 e2_flip = -1;
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

          Real cx00 = old_coords[v0_old_face*dim + 0];
          Real cy00 = old_coords[v0_old_face*dim + 1];
          //Real cz00 = old_coords[v0_old_face*dim + 2];
          Real cx30 = old_coords[v1*dim + 0];
          Real cy30 = old_coords[v1*dim + 1];
          //Real cz30 = old_coords[v1*dim + 2];
          Real cx03 = old_coords[v2*dim + 0];
          Real cy03 = old_coords[v2*dim + 1];
          //Real cz03 = old_coords[v2*dim + 2];

          auto pts_per_edge = n_edge_pts;
          Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
          Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
          //Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
          Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
          Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e0_flip > 0) {
            auto tempx = cx10;
            auto tempy = cy10;
            //auto tempz = cz10;
            cx10 = cx20;
            cy10 = cy20;
            //cz10 = cz20;
            cx20 = tempx;
            cy20 = tempy;
            //cz20 = tempz;
          }

          Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
          Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
          //Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
          Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e1_flip > 0) {
            auto tempx = cx21;
            auto tempy = cy21;
            //auto tempz = cz21;
            cx21 = cx12;
            cy21 = cy12;
            //cz21 = cz12;
            cx12 = tempx;
            cy12 = tempy;
            //cz12 = tempz;
          }

          Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
          Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
          //Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
          Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e2_flip > 0) {
            auto tempx = cx02;
            auto tempy = cy02;
            //auto tempz = cz02;
            cx02 = cx01;
            cy02 = cy01;
            //cz02 = cz01;
            cx01 = tempx;
            cy01 = tempy;
            //cz01 = tempz;
          }

          Real cx11 = old_faceCtrlPts[old_face*dim + 0];
          Real cy11 = old_faceCtrlPts[old_face*dim + 1];
          //Real cz11 = old_faceCtrlPts[old_face*dim + 2];

          //get the interp points
          auto p1_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
                      cx10*B10_cube(nodePts[0], nodePts[1]) +
                      cx20*B20_cube(nodePts[0], nodePts[1]) +
                      cx30*B30_cube(nodePts[0], nodePts[1]) +
                      cx21*B21_cube(nodePts[0], nodePts[1]) +
                      cx12*B12_cube(nodePts[0], nodePts[1]) +
                      cx03*B03_cube(nodePts[0], nodePts[1]) +
                      cx02*B02_cube(nodePts[0], nodePts[1]) +
                      cx01*B01_cube(nodePts[0], nodePts[1]) +
                      cx11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
                      cy10*B10_cube(nodePts[0], nodePts[1]) +
                      cy20*B20_cube(nodePts[0], nodePts[1]) +
                      cy30*B30_cube(nodePts[0], nodePts[1]) +
                      cy21*B21_cube(nodePts[0], nodePts[1]) +
                      cy12*B12_cube(nodePts[0], nodePts[1]) +
                      cy03*B03_cube(nodePts[0], nodePts[1]) +
                      cy02*B02_cube(nodePts[0], nodePts[1]) +
                      cy01*B01_cube(nodePts[0], nodePts[1]) +
                      cy11*B11_cube(nodePts[0], nodePts[1]);
          auto p2_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
                      cx10*B10_cube(nodePts[2], nodePts[3]) +
                      cx20*B20_cube(nodePts[2], nodePts[3]) +
                      cx30*B30_cube(nodePts[2], nodePts[3]) +
                      cx21*B21_cube(nodePts[2], nodePts[3]) +
                      cx12*B12_cube(nodePts[2], nodePts[3]) +
                      cx03*B03_cube(nodePts[2], nodePts[3]) +
                      cx02*B02_cube(nodePts[2], nodePts[3]) +
                      cx01*B01_cube(nodePts[2], nodePts[3]) +
                      cx11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
                      cy10*B10_cube(nodePts[2], nodePts[3]) +
                      cy20*B20_cube(nodePts[2], nodePts[3]) +
                      cy30*B30_cube(nodePts[2], nodePts[3]) +
                      cy21*B21_cube(nodePts[2], nodePts[3]) +
                      cy12*B12_cube(nodePts[2], nodePts[3]) +
                      cy03*B03_cube(nodePts[2], nodePts[3]) +
                      cy02*B02_cube(nodePts[2], nodePts[3]) +
                      cy01*B01_cube(nodePts[2], nodePts[3]) +
                      cy11*B11_cube(nodePts[2], nodePts[3]);

          //use these as interp pts to find ctrl pts in new mesh
          {
            Real cx0 = old_coords[old_vert_noKey*dim + 0];
            Real cy0 = old_coords[old_vert_noKey*dim + 1];

            Matrix<2,1> cx({cx0, vert_ctrlPts[mid_vert*1*dim + 0]});
            Matrix<2,1> cy({cy0, vert_ctrlPts[mid_vert*1*dim + 1]});

            Matrix<2,1> fx({p1_x, p2_x});
            Matrix<2,1> fy({p1_y, p2_y});

            Matrix<2,2> M1_inv({B1_cube(xi_1_cube()), B2_cube(xi_1_cube()), B1_cube(xi_2_cube()),
                B2_cube(xi_2_cube())});
            Matrix<2,2> M2({B0_cube(xi_1_cube()), B3_cube(xi_1_cube()), B0_cube(xi_2_cube()),
                B3_cube(xi_2_cube())});

            auto M1 = invert(M1_inv);
            auto Cx = M1*fx - M1*M2*cx;
            auto Cy = M1*fy - M1*M2*cy;
            edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
          }
        }
        //edge_curveXi[new_e2*2 + 0] = 0.0;
        //edge_curveXi[new_e2*2 + 1] = 1.0;

      }
      atomic_fetch_add(&count_key[0], 1);
    }
  };
  parallel_for(nold_edge, std::move(create_crv_edges));

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim);
  new_mesh->add_tag<Real>(0, "bezier_pts", dim);
  //new_mesh->add_tag<Real>(1, "edge_ends_xi", 2, edge_curveXi);
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));
  new_mesh->set_tag_for_ctrlPts(0, Reals(vert_ctrlPts));

  return;
}

/*
//TODO transfer face pts
void create_curved_faces(Mesh *mesh, Mesh *new_mesh, LOs old2new, LOs prods2new,
                         LOs keys2prods, LOs keys2midverts,
                         LOs old_verts2new_verts) {
  OMEGA_H_TIME_FUNCTION;
  auto nold_edge = old2new.size();
  auto old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto old_fe2e = mesh->get_adj(2, 1).ab2b;
  auto old_ef2f = mesh->ask_up(1, 2).ab2b;
  auto old_e2ef = mesh->ask_up(1, 2).a2ab;
  auto old_fv2v = mesh->ask_down(2, 0).ab2b;
  auto old_coords = mesh->coords();
  auto old_edgeCtrlPts = mesh->get_ctrlPts(1);
  auto old_faceCtrlPts = mesh->get_ctrlPts(2);
  auto order = mesh->get_max_order();
  auto dim = mesh->dim();
  auto n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto n_face_pts = mesh->n_internal_ctrlPts(2);

  auto new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto nnew_face = new_mesh->nfaces();
  auto nnew_verts = new_mesh->nverts();

  Write<Real> face_ctrlPts(nnew_face*n_face_pts*dim);

  auto new_verts2old_verts = invert_map_by_atomics(old_verts2new_verts,
                                                   nnew_verts);

  Write<LO> count_key(1, 0);
  OMEGA_H_CHECK(order == 3);

  auto create_crv_edges = OMEGA_H_LAMBDA (LO old_face) {
    LO v0_old = old_ev2v[old_face*2 + 0];
    LO v1_old = old_ev2v[old_face*2 + 1];
    LO v2_old = old_ev2v[old_face*2 + 2];

    if (old2new[old_face] != -1) {
      // map as is from old id to new id
      LO new_face = old2new[old_face];
      printf("old edge %d same as new edge %d\n", old_face, new_face);
      for (I8 d = 0; d < dim; ++d) {
        face_ctrlPts[new_face*n_face_pts*dim + d] =
          old_faceCtrlPts[old_face*n_face_pts*dim + d];
        face_ctrlPts[new_face*n_face_pts*dim + (n_face_pts-1)*dim + d] =
          old_faceCtrlPts[old_face*n_face_pts*dim + (n_face_pts-1)*dim + d];
      }
    }
    else {
      LO key_id = count_key[0];
      LO mid_vert = keys2midverts[key_id];
      LO start = keys2prods[key_id];
      LO end = keys2prods[key_id + 1] - 1;
      OMEGA_H_CHECK((end-start) == 1);
      LO new_e0 = prods2new[start];
      LO new_e1 = prods2new[start+1];
      LO new_e2 = prods2new[start+2];

      LO v1_new_e0 = new_ev2v[new_e0*2 + 1];
      LO v0_new_e1 = new_ev2v[new_e1*2 + 0];
      OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));
      //LO v0_new_e0 = new_ev2v[new_e0*2 + 0];
      //LO v1_new_e1 = new_ev2v[new_e1*2 + 1];
      //OMEGA_H_CHECK((v0_new_e0 == v0_old) && (v1_new_e1 == v1_old));
      //topologically yes but diff mesh objects

      LO v0_new_e2 = new_ev2v[new_e2*2 + 0];
      LO v1_new_e2 = new_ev2v[new_e2*2 + 1];
      OMEGA_H_CHECK(v1_new_e2 == mid_vert);
      // swap unneeded, order checked

      // check flip
      //auto flip_new_e0 = -1;
      //auto flip_new_e1 = -1;

      //ctrl pts for e0
      {
        //for 2d mesh for now, order=3
        //query old ctrl pts
        Real cx0 = old_coords[v0_old*dim + 0];
        Real cy0 = old_coords[v0_old*dim + 1];
        Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 0];
        Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim +
          (n_edge_pts-1)*dim + 1];
        Real cx3 = old_coords[v1_old*dim + 0];
        Real cy3 = old_coords[v1_old*dim + 1];

        Real new_xi_3 = 0.5;
        Real new_cx3 = cx0*B0_cube(new_xi_3) + cx1*B1_cube(new_xi_3) +
                       cx2*B2_cube(new_xi_3) + cx3*B3_cube(new_xi_3);
        Real new_cy3 = cy0*B0_cube(new_xi_3) + cy1*B1_cube(new_xi_3) +
                       cy2*B2_cube(new_xi_3) + cy3*B3_cube(new_xi_3);
        printf(" new mid interp pt %f, %f \n", new_cx3, new_cy3);
        vert_ctrlPts[mid_vert*1*dim + 0] = new_cx3;
        vert_ctrlPts[mid_vert*1*dim + 1] = new_cy3;
        //vert_ctrlPts[mid_vert*1*dim + 2] = new_cz3;

        Real old_xi_2 = xi_2_cube();
        Real new_xi_2 = old_xi_2/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);
        printf("edgo 0 p2 %f, %f , xi2 %f\n", new_px2, new_py2, new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = old_xi_1/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);
        printf("edgo 0 p1 %f, %f \n", new_px1, new_py1);

        Matrix<2,1> c0({cx0, cy0});
        Matrix<2,1> c3({new_cx3, new_cy3});

        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        auto M1 = invert(M1_inv);
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        auto Cx = M1*fx - M1*M2*c0;
        auto Cy = M1*fy - M1*M2*c3;
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        //edge_curveXi[new_e0*2 + 0] = 0.0;
        //edge_curveXi[new_e0*2 + 1] = 1.0;
      }
      //ctrl pts for e1
      {
        //for 2d mesh for now, order=3
        //query old ctrl pts
        Real cx0 = old_coords[v0_old*dim + 0];
        Real cy0 = old_coords[v0_old*dim + 1];
        Real cx1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 0];
        Real cy1 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + 1];
        Real cx2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 0];
        Real cy2 = old_edgeCtrlPts[old_edge*n_edge_pts*dim + (n_edge_pts-1)*dim + 1];
        Real cx3 = old_coords[v1_old*dim + 0];
        Real cy3 = old_coords[v1_old*dim + 1];

        Real old_xi_2 = xi_2_cube();
        Real new_xi_2 = 0.5 + old_xi_2/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);

        Real old_xi_1 = xi_1_cube();
        Real new_xi_1 = 0.5 + old_xi_1/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);

        Matrix<2,1> c0({vert_ctrlPts[mid_vert*1*dim + 0],
                        vert_ctrlPts[mid_vert*1*dim + 1]});
        Matrix<2,1> c3({cx3, cy3});

        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({B1_cube(old_xi_1), B2_cube(old_xi_1), B1_cube(old_xi_2),
                            B2_cube(old_xi_2)});
        auto M1 = invert(M1_inv);
        Matrix<2,2> M2({B0_cube(old_xi_1), B3_cube(old_xi_1), B0_cube(old_xi_2),
                        B3_cube(old_xi_2)});

        auto Cx = M1*fx - M1*M2*c0;
        auto Cy = M1*fy - M1*M2*c3;
        edge_ctrlPts[new_e1*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e1*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
        //edge_curveXi[new_e1*2 + 0] = 0.0;
        //edge_curveXi[new_e1*2 + 1] = 1.0;
      }
      //ctrl pts for e2
      {
        auto ab2b = new_verts2old_verts.ab2b;
        auto a2ab = new_verts2old_verts.a2ab;
        LO old_vert_noKey = ab2b[a2ab[v0_new_e2]];

        LO old_face = -1;
        for (LO index = old_e2ef[old_edge]; index < old_e2ef[old_edge + 1];
             ++index) {
          LO adj_face = old_ef2f[index];
          for (LO vert = 0; vert < 3; ++vert) {
            LO vert_old_face = old_fv2v[adj_face*3 + vert];
            if (vert_old_face == old_vert_noKey) {
              old_face = adj_face;
              break;
            }
          }
          if (old_face > 0) {
            break;
          }
        }

        // for that old-face, query points at calculated xi values
        printf("For old edge %d, found old face %d\n", old_edge, old_face);
        {
          auto v0_old_face = old_fv2v[old_face*3];
          auto v1 = old_fv2v[old_face*3 + 1];
          auto v2 = old_fv2v[old_face*3 + 2];

          //TODO check alignment of new noKey edge
          //from the example, for now p1 is nearer to noKey vert
          auto nodePts = xi_values_cubic(old_vert_noKey, v0_old_face, v1, v2);

          auto old_face_e0 = old_fe2e[old_face*3];
          auto old_face_e1 = old_fe2e[old_face*3 + 1];
          auto old_face_e2 = old_fe2e[old_face*3 + 2];

          I8 e0_flip = -1;
          I8 e1_flip = -1;
          I8 e2_flip = -1;
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

          Real cx00 = old_coords[v0_old_face*dim + 0];
          Real cy00 = old_coords[v0_old_face*dim + 1];
          //Real cz00 = old_coords[v0_old_face*dim + 2];
          Real cx30 = old_coords[v1*dim + 0];
          Real cy30 = old_coords[v1*dim + 1];
          //Real cz30 = old_coords[v1*dim + 2];
          Real cx03 = old_coords[v2*dim + 0];
          Real cy03 = old_coords[v2*dim + 1];
          //Real cz03 = old_coords[v2*dim + 2];

          auto pts_per_edge = n_edge_pts;
          Real cx10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 0];
          Real cy10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 1];
          //Real cz10 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + 2];
          Real cx20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];//2 pts per edge
          Real cy20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz20 = old_edgeCtrlPts[old_face_e0*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e0_flip > 0) {
            auto tempx = cx10;
            auto tempy = cy10;
            //auto tempz = cz10;
            cx10 = cx20;
            cy10 = cy20;
            //cz10 = cz20;
            cx20 = tempx;
            cy20 = tempy;
            //cz20 = tempz;
          }

          Real cx21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 0];
          Real cy21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 1];
          //Real cz21 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + 2];
          Real cx12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz12 = old_edgeCtrlPts[old_face_e1*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e1_flip > 0) {
            auto tempx = cx21;
            auto tempy = cy21;
            //auto tempz = cz21;
            cx21 = cx12;
            cy21 = cy12;
            //cz21 = cz12;
            cx12 = tempx;
            cy12 = tempy;
            //cz12 = tempz;
          }

          Real cx02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 0];
          Real cy02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 1];
          //Real cz02 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + 2];
          Real cx01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 0];
          Real cy01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 1];
          //Real cz01 = old_edgeCtrlPts[old_face_e2*pts_per_edge*dim + (pts_per_edge-1)*dim + 2];
          if (e2_flip > 0) {
            auto tempx = cx02;
            auto tempy = cy02;
            //auto tempz = cz02;
            cx02 = cx01;
            cy02 = cy01;
            //cz02 = cz01;
            cx01 = tempx;
            cy01 = tempy;
            //cz01 = tempz;
          }

          Real cx11 = old_faceCtrlPts[old_face*dim + 0];
          Real cy11 = old_faceCtrlPts[old_face*dim + 1];
          //Real cz11 = old_faceCtrlPts[old_face*dim + 2];

          //get the interp points
          auto p1_x = cx00*B00_cube(nodePts[0], nodePts[1]) +
                      cx10*B10_cube(nodePts[0], nodePts[1]) +
                      cx20*B20_cube(nodePts[0], nodePts[1]) +
                      cx30*B30_cube(nodePts[0], nodePts[1]) +
                      cx21*B21_cube(nodePts[0], nodePts[1]) +
                      cx12*B12_cube(nodePts[0], nodePts[1]) +
                      cx03*B03_cube(nodePts[0], nodePts[1]) +
                      cx02*B02_cube(nodePts[0], nodePts[1]) +
                      cx01*B01_cube(nodePts[0], nodePts[1]) +
                      cx11*B11_cube(nodePts[0], nodePts[1]);
          auto p1_y = cy00*B00_cube(nodePts[0], nodePts[1]) +
                      cy10*B10_cube(nodePts[0], nodePts[1]) +
                      cy20*B20_cube(nodePts[0], nodePts[1]) +
                      cy30*B30_cube(nodePts[0], nodePts[1]) +
                      cy21*B21_cube(nodePts[0], nodePts[1]) +
                      cy12*B12_cube(nodePts[0], nodePts[1]) +
                      cy03*B03_cube(nodePts[0], nodePts[1]) +
                      cy02*B02_cube(nodePts[0], nodePts[1]) +
                      cy01*B01_cube(nodePts[0], nodePts[1]) +
                      cy11*B11_cube(nodePts[0], nodePts[1]);
          auto p2_x = cx00*B00_cube(nodePts[2], nodePts[3]) +
                      cx10*B10_cube(nodePts[2], nodePts[3]) +
                      cx20*B20_cube(nodePts[2], nodePts[3]) +
                      cx30*B30_cube(nodePts[2], nodePts[3]) +
                      cx21*B21_cube(nodePts[2], nodePts[3]) +
                      cx12*B12_cube(nodePts[2], nodePts[3]) +
                      cx03*B03_cube(nodePts[2], nodePts[3]) +
                      cx02*B02_cube(nodePts[2], nodePts[3]) +
                      cx01*B01_cube(nodePts[2], nodePts[3]) +
                      cx11*B11_cube(nodePts[2], nodePts[3]);
          auto p2_y = cy00*B00_cube(nodePts[2], nodePts[3]) +
                      cy10*B10_cube(nodePts[2], nodePts[3]) +
                      cy20*B20_cube(nodePts[2], nodePts[3]) +
                      cy30*B30_cube(nodePts[2], nodePts[3]) +
                      cy21*B21_cube(nodePts[2], nodePts[3]) +
                      cy12*B12_cube(nodePts[2], nodePts[3]) +
                      cy03*B03_cube(nodePts[2], nodePts[3]) +
                      cy02*B02_cube(nodePts[2], nodePts[3]) +
                      cy01*B01_cube(nodePts[2], nodePts[3]) +
                      cy11*B11_cube(nodePts[2], nodePts[3]);

          //use these as interp pts to find ctrl pts in new mesh
          {
            Real cx0 = old_coords[old_vert_noKey*dim + 0];
            Real cy0 = old_coords[old_vert_noKey*dim + 1];

            Matrix<2,1> c0({cx0, cy0});
            Matrix<2,1> c3({vert_ctrlPts[mid_vert*1*dim + 0],
                            vert_ctrlPts[mid_vert*1*dim + 1]});

            Matrix<2,1> fx({p1_x, p2_x});
            Matrix<2,1> fy({p1_y, p2_y});

            Matrix<2,2> M1_inv({B1_cube(xi_1_cube()), B2_cube(xi_1_cube()), B1_cube(xi_2_cube()),
                B2_cube(xi_2_cube())});
            auto M1 = invert(M1_inv);
            Matrix<2,2> M2({B0_cube(xi_1_cube()), B3_cube(xi_1_cube()), B0_cube(xi_2_cube()),
                B3_cube(xi_2_cube())});

            auto Cx = M1*fx - M1*M2*c0;
            auto Cy = M1*fy - M1*M2*c3;
            edge_ctrlPts[new_e2*n_edge_pts*dim + 0] = Cx(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + 1] = Cy(0,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
            edge_ctrlPts[new_e2*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);
          }
        }
        //edge_curveXi[new_e2*2 + 0] = 0.0;
        //edge_curveXi[new_e2*2 + 1] = 1.0;

      }
      atomic_fetch_add(&count_key[0], 1);
    }
  };
  parallel_for(nold_edge, std::move(create_crv_edges));

  new_mesh->add_tag<Real>(1, "bezier_pts", n_edge_pts*dim);
  new_mesh->add_tag<Real>(0, "bezier_pts", dim);
  //new_mesh->add_tag<Real>(1, "edge_ends_xi", 2, edge_curveXi);
  new_mesh->set_tag_for_ctrlPts(1, Reals(edge_ctrlPts));
  new_mesh->set_tag_for_ctrlPts(0, Reals(vert_ctrlPts));

  return;
}
*/

#define OMEGA_H_INST(T)
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h
