#include "Omega_h_curve_coarsen.hpp"

#include "Omega_h_curve_validity_3d.hpp"

namespace Omega_h {

void check_validity_all_tet(Mesh *new_mesh) {

  auto const new_rv2v = new_mesh->ask_down(3, 0).ab2b;
  auto const new_re2e = new_mesh->ask_down(3, 1).ab2b;
  auto const new_rf2f = new_mesh->ask_down(3, 2).ab2b;
  auto const new_ev2v = new_mesh->ask_down(1, 0).ab2b;
  auto const new_e2r = new_mesh->ask_up(1, 3);
  auto const new_e2er = new_e2r.a2ab;
  auto const new_er2r = new_e2r.ab2b;
  auto const nnew_tet = new_mesh->nregions();

  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const faceCtrlPts = new_mesh->get_ctrlPts(2);

  printf("mesh size %d tets\n", nnew_tet);
  auto check_tet = OMEGA_H_LAMBDA(LO i) {
    LO has_invalid_tet = -1;
    Few<Real, 60> tet_pts = collect_tet_pts(3, i, new_ev2v, new_rv2v, vertCtrlPts
        , edgeCtrlPts, faceCtrlPts, new_re2e, new_rf2f);

    Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

    has_invalid_tet = checkMinJacDet_3d(nodes_det, 3);
    if (has_invalid_tet > 0) printf("tet %d invalid code %d\n", i, has_invalid_tet);
  };
  parallel_for(nnew_tet, std::move(check_tet), "check_tet");
  return;
}

void correct_curved_edges(Mesh *new_mesh) {

  auto const edge_crv2bdry_dim = new_mesh->get_array<I8>(1, "edge_crv2bdry_dim");
  auto const new_rv2v = new_mesh->ask_down(3, 0).ab2b;
  auto const new_re2e = new_mesh->ask_down(3, 1).ab2b;
  auto const new_rf2f = new_mesh->ask_down(3, 2).ab2b;
  auto const new_ev2v = new_mesh->ask_down(1, 0).ab2b;
  auto const new_e2r = new_mesh->ask_up(1, 3);
  auto const new_e2er = new_e2r.a2ab;
  auto const new_er2r = new_e2r.ab2b;
  auto const new_e2f = new_mesh->ask_up(1, 2);
  auto const new_e2ef = new_e2f.a2ab;
  auto const new_ef2f = new_e2f.ab2b;
  auto const nnew_edge = new_mesh->nedges();
  auto const nnew_face = new_mesh->nfaces();
  auto const n_edge_pts = new_mesh->n_internal_ctrlPts(1);

  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const faceCtrlPts = new_mesh->get_ctrlPts(2);

  Write<I8> face_vis (nnew_face, -1);
  //pfor over all edges using curved to bdry algo
  //if it was curved before, check validity of all its adjacent tets
  //if any adjacent tet is invalid, print out error
  auto edge_correct = OMEGA_H_LAMBDA(LO i) {
    LO has_invalid_tet = -1;
    if (edge_crv2bdry_dim[i] >= 1) {
      for (LO er = new_e2er[i]; er < new_e2er[i+1]; ++er) {
        if (has_invalid_tet > 0) break;
        LO adj_tet = new_er2r[er];
        Few<Real, 60> tet_pts = collect_tet_pts(3, adj_tet, new_ev2v, new_rv2v, vertCtrlPts
            , edgeCtrlPts, faceCtrlPts, new_re2e, new_rf2f);

        Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

        has_invalid_tet = checkMinJacDet_3d(nodes_det, n_edge_pts+1);
        if (has_invalid_tet > 0) {
          printf(" invalid tet %d for edge %d curved to g_dim %d with code %d\n"
              , adj_tet, i, edge_crv2bdry_dim[i], has_invalid_tet);
          for (LO k=0; k<4; ++k) face_vis[new_rf2f[adj_tet*4 + k]] = 1;
        }
      }
      if (has_invalid_tet > 0) {
        for (LO ef = new_e2ef[i]; ef < new_e2ef[i+1]; ++ef) {
          //LO adj_f = new_ef2f[ef];
          //face_vis[adj_f] = 1;
        }
      }
    }
  };
  parallel_for(nnew_edge, std::move(edge_correct), "edge_correct");
  new_mesh->set_tag<I8>(2, "face_dualCone", Read<I8>(face_vis));
  return;
}

void check_validity_new_curved_edges(Mesh *new_mesh) {

  auto const new_rv2v = new_mesh->ask_down(3, 0).ab2b;
  auto const new_re2e = new_mesh->ask_down(3, 1).ab2b;
  auto const new_rf2f = new_mesh->ask_down(3, 2).ab2b;
  auto const new_ev2v = new_mesh->ask_down(1, 0).ab2b;
  auto const new_e2r = new_mesh->ask_up(1, 3);
  auto const new_e2er = new_e2r.a2ab;
  auto const new_er2r = new_e2r.ab2b;
  auto const nnew_edge = new_mesh->nedges();
  auto const n_edge_pts = new_mesh->n_internal_ctrlPts(1);

  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const faceCtrlPts = new_mesh->get_ctrlPts(2);
  auto const edge_dualCone = new_mesh->get_array<I8>(1, "edge_dualCone");

  //pfor over all edges of dual cone cav
  //check validity of all its adjacent tets
  auto edge_correct = OMEGA_H_LAMBDA(LO i) {
    LO has_invalid_tet = -1;
    if (edge_dualCone[i] == 1) {
      for (LO er = new_e2er[i]; er < new_e2er[i+1]; ++er) {
        LO adj_tet = new_er2r[er];
        Few<Real, 60> tet_pts = collect_tet_pts(3, adj_tet, new_ev2v, new_rv2v, vertCtrlPts
            , edgeCtrlPts, faceCtrlPts, new_re2e, new_rf2f);

        Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

        has_invalid_tet = checkMinJacDet_3d(nodes_det, n_edge_pts+1);
        if (has_invalid_tet > 0) {
          printf("dualCone edge %d in interior creates invalid tet with code %d\n",
              i, has_invalid_tet);
          break;
        }
      }
    }
  };
  parallel_for(nnew_edge, std::move(edge_correct), "edge_correct");
  return;
}

} // namespace Omega_h
