#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

LOs checkValidity_2d(Mesh *mesh, LOs new_tris, Int const mesh_dim) {
  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto const n_edge_pts = mesh->n_internal_ctrlPts(1);
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order == 3);

  Write<LO> is_invalid(new_tris.size(), -1);
  //LO const ntri_pts = 10;

  auto check_validity = OMEGA_H_LAMBDA (LO n) {
    //TODO change tri_pts size to 30 for 3d
    Few<Real, 20> tri_pts;//ntri_pts*dim=20
    auto tri = new_tris[n];

    //query the tri's down verts's ctrl pts and store
    for (LO j = 0; j < 3; ++j) {
      if (mesh_dim == 2) {
        auto p = get_vector<2>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
      else {
        OMEGA_H_CHECK (mesh_dim == 3);
        auto p = get_vector<3>(vertCtrlPts, fv2v[tri*3 + j]);
        for (LO k = 0; k < mesh_dim; ++k) {
          tri_pts[j*mesh_dim + k] = p[k];
        }
      }
    }

    //query the tri's down edge's ctrl pts and store

    auto v0 = fv2v[tri*3 + 0];
    auto v1 = fv2v[tri*3 + 1];
    auto v2 = fv2v[tri*3 + 2];
    auto e0 = fe2e[tri*3 + 0];
    auto e1 = fe2e[tri*3 + 1];
    auto e2 = fe2e[tri*3 + 2];
    auto e0v0 = ev2v[e0*2 + 0];
    auto e0v1 = ev2v[e0*2 + 1];
    auto e1v0 = ev2v[e1*2 + 0];
    auto e1v1 = ev2v[e1*2 + 1];
    auto e2v0 = ev2v[e2*2 + 0];
    auto e2v1 = ev2v[e2*2 + 1];
    auto flip = vector_3(-1, -1, -1);
    if ((e0v0 == v1) && (e0v1 == v0)) {
      flip[0] = 1;
    }
    else {
      OMEGA_H_CHECK((e0v0 == v0) && (e0v1 == v1));
    }
    if ((e1v0 == v2) && (e1v1 == v1)) {
      flip[1] = 1;
    }
    else {
      OMEGA_H_CHECK((e1v0 == v1) && (e1v1 == v2));
    }
    if ((e2v0 == v0) && (e2v1 == v2)) {
      flip[2] = 1;
    }

    for (LO j = 0; j < 3; ++j) {
      LO index = 3;
      if (flip[j] == -1) {
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
        }
      }
      else {
        //for flipped edges
        OMEGA_H_CHECK (flip[j] == 1);
        for (I8 d = 0; d < mesh_dim; ++d) {
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + mesh_dim + d];
          tri_pts[index*mesh_dim + j*n_edge_pts*mesh_dim + mesh_dim + d] =
            edgeCtrlPts[fe2e[tri*3 + j]*n_edge_pts*mesh_dim + d];
        }
      }
    }

    //query the face's ctrl pt and store
    for (I8 d = 0; d < mesh_dim; ++d) {
      LO index = 9;
      tri_pts[index*mesh_dim + d] = faceCtrlPts[tri*mesh_dim + d];
    }

    //TODO change to template for mesh_dim
    auto nodes_det = getTriJacDetNodes<15, 2>(order, tri_pts);

    is_invalid[n] = checkMinJacDet<15>(nodes_det, order);
    if (is_invalid[n] > 0) printf("invalid tri %d\n", tri);
  };
  parallel_for(new_tris.size(), std::move(check_validity));

  return LOs(is_invalid);
}

} //namespace
