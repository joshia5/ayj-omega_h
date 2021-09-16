#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_atomics.hpp"

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

    Real cx0 = coords[v0*dim + 0];
    Real cy0 = coords[v0*dim + 1];
    Real cz0 = coords[v0*dim + 2];
    Real x1 = interpPts[i*dim + 0];
    Real y1 = interpPts[i*dim + 1];
    Real z1 = interpPts[i*dim + 2];
    Real cx2 = coords[v1*dim + 0];
    Real cy2 = coords[v1*dim + 1];
    Real cz2 = coords[v1*dim + 2];

    auto cx1 = (x1 - B0_quad(xi_1)*cx0 - B2_quad(xi_1)*cx2)/B1_quad(xi_1);
    auto cy1 = (y1 - B0_quad(xi_1)*cy0 - B2_quad(xi_1)*cy2)/B1_quad(xi_1);
    auto cz1 = (z1 - B0_quad(xi_1)*cz0 - B2_quad(xi_1)*cz2)/B1_quad(xi_1);
    new_pts[i*dim + 0] = cx1;
    new_pts[i*dim + 1] = cy1;
    new_pts[i*dim + 2] = cz1;
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

void transfer_bezier_edges(Mesh *mesh, Mesh *new_mesh,
    LOs old2new, LOs prods2new, LOs keys2prods,
    LOs keys2midverts, LOs keys2edges) {
  
  LO count_key = 0;
  auto nold_edge = old2new.size();
  auto old_ev2v = mesh->get_adj(1, 0).ab2b;
  auto new_ev2v = new_mesh->get_adj(1, 0).ab2b;
  auto old_coords = mesh->coords();
  auto new_coords = new_mesh->coords();
  auto n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto order = mesh->get_max_order;
  auto dim = mesh->dim();
  auto old_ctrlPts = mesh->get_ctrlPts(1);
  auto nnew_edge = new_ev2v.size()/2;
  Write<Real> edge_ctrlPts(nnew_edge*n_edge_pts*dim);

  auto transfer_edges = OMEGA_H_LAMBDA (LO old_edge) {
    v0_old = old_ev2v[old_edge*2 + 0];
    v1_old = old_ev2v[old_edge*2 + 1];

    if (old2new[i] != -1) {
      new_edge = old2new[i];
      // map as is from old id to new id
    }
    else {
      key_id = count_key;
      mid_vert = keys2midverts[key_id];
      start = prods2new[i];
      end = prods2new[i+1] - 1;
      OMEGA_H_CHECK((end-start) == 2);
      new_e0 = prods2new[start];
      new_e1 = prods2new[start+1];
      new_e2 = prods2new[start+2];

      v0_new_e0 = new_ev2v[new_e0*2 + 0];
      v1_new_e0 = new_ev2v[new_e0*2 + 1];
      v0_new_e1 = new_ev2v[new_e1*2 + 0];
      v1_new_e1 = new_ev2v[new_e1*2 + 1];
      OMEGA_H_CHECK((v1_new_e0 == mid_vert) && (v0_new_e1 == mid_vert));
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
        Real cx0 = new_coords[v0_old*dim + 0];
        Real cy0 = new_coords[v0_old*dim + 1];
        Real cx1 = old_ctrlPts[i*n_edge_pts*dim + 0];
        Real cy1 = old_ctrlPts[i*n_edge_pts*dim + 1];
        Real cx2 = old_ctrlPts[i*n_edge_pts*dim + (n_edge_pts-1)*dim + 0];
        Real cy2 = old_ctrlPts[i*n_edge_pts*dim + (n_edge_pts-1)*dim + 1];
        Real cx3 = new_coords[v1_old*dim + 0];
        Real cy3 = new_coords[v1_old*dim + 1];

        Real new_xi_3 = 0.5;
        Real new_cx3 = cx0*B0_cube(new_xi_3) + cx1*B1_cube(new_xi_3) +
                       cx2*B2_cube(new_xi_3) + cx3*B3_cube(new_xi_3);
        Real new_cy3 = cy0*B0_cube(new_xi_3) + cy1*B1_cube(new_xi_3) +
                       cy2*B2_cube(new_xi_3) + cy3*B3_cube(new_xi_3);

        Real new_xi_2 = xi_2_cube()/2.0;
        Real new_px2 = cx0*B0_cube(new_xi_2) + cx1*B1_cube(new_xi_2) +
                       cx2*B2_cube(new_xi_2) + cx3*B3_cube(new_xi_2);
        Real new_py2 = cy0*B0_cube(new_xi_2) + cy1*B1_cube(new_xi_2) +
                       cy2*B2_cube(new_xi_2) + cy3*B3_cube(new_xi_2);

        Real new_xi_1 = xi_1_cube()/2.0;
        Real new_px1 = cx0*B0_cube(new_xi_1) + cx1*B1_cube(new_xi_1) +
                       cx2*B2_cube(new_xi_1) + cx3*B3_cube(new_xi_1);
        Real new_py1 = cy0*B0_cube(new_xi_1) + cy1*B1_cube(new_xi_1) +
                       cy2*B2_cube(new_xi_1) + cy3*B3_cube(new_xi_1);

        Matrix<2,1> c0({cx0, cy0});
        Matrix<2,1> c3({new_cx3, new_cy3});

        Matrix<2,1> fx({new_px1, new_px2});
        Matrix<2,1> fy({new_py1, new_py2});

        Matrix<2,2> M1_inv({B1_cube(new_xi_1), B2_cube(new_xi_1), B1_cube(new_xi_2),
                            B2_cube(new_xi_2)});
        auto M1 = invert(M1_inv);
        Matrix<2,2> M2({B0_cube(new_xi_1), B3_cube(new_xi_1), B0_cube(new_xi_2),
                        B3_cube(new_xi_2)});

        auto Cx = M1*fx - M1*M2*p0;
        auto Cy = M1*fy - M1*M2*p1;
        edge_ctrlPts[new_e0*n_edge_pts*dim + 0] = Cx(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + 1] = Cy(0,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 0] = Cx(1,0);
        edge_ctrlPts[new_e0*n_edge_pts*dim + (n_edge_pts-1)*dim + 1] = Cy(1,0);

      }
      //ctrl pts for e1
      {
      }
      //ctrl pts for e2
      {
        inquire old face
      }

      atomic_fetch_add(&count_key, 1);
    }
  };
  parallel_for(nold_edge, std::move(transfer_edges));

  return;
}

#define OMEGA_H_INST(T)
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

} // namespace Omega_h
