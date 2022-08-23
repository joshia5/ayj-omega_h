#include "Omega_h_curve_coarsen.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_file.hpp"

#include "Omega_h_curve_validity_3d.hpp"

namespace Omega_h {

void check_validity_all_tet(Mesh *new_mesh) {

  auto const new_ev2v = new_mesh->ask_down(1, 0).ab2b;
  auto const new_rv2v = new_mesh->ask_down(3, 0).ab2b;
  auto const new_re2e = new_mesh->ask_down(3, 1).ab2b;
  auto const new_rf2f = new_mesh->ask_down(3, 2).ab2b;
  auto const new_e2r = new_mesh->ask_up(1, 3);
  auto const new_e2er = new_e2r.a2ab;
  auto const new_er2r = new_e2r.ab2b;
  auto const nnew_tet = new_mesh->nregions();
  auto const nfaces = new_mesh->nfaces();

  auto const vertCtrlPts = new_mesh->get_ctrlPts(0);
  auto const edgeCtrlPts = new_mesh->get_ctrlPts(1);
  auto const faceCtrlPts = new_mesh->get_ctrlPts(2);

  Write<LO> invalid_tet(nnew_tet, -1);

  printf("mesh size %d tets\n", nnew_tet);
  auto check_tet = OMEGA_H_LAMBDA(LO i) {
    LO is_invalid = -1;
    Few<Real, 60> tet_pts = collect_tet_pts(3, i, new_ev2v, new_rv2v, vertCtrlPts
        , edgeCtrlPts, faceCtrlPts, new_re2e, new_rf2f);

    Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

    is_invalid = checkMinJacDet_3d(nodes_det, 3, -1);
    invalid_tet[i] = is_invalid;
    if (is_invalid > 0) {
      printf("tet %d invalid code %d\n", i, is_invalid);
    }
  };
  parallel_for(nnew_tet, std::move(check_tet), "check_tet");

  new_mesh->add_tag<LO>(3, "invalidity", 1, Read<LO>(invalid_tet));
  HostWrite<LO> face_crvVis(nfaces);
  for (LO i = 0; i < nfaces; ++i) {
    face_crvVis[i] = -1;
  }
  if (new_mesh->has_tag(2, "face_crvVis")) {
    new_mesh->remove_tag(2, "face_crvVis");
  }
  new_mesh->add_tag<LO>(2, "face_crvVis", 1);
  auto tet_invalid_h = HostRead<LO>(Read<LO>(invalid_tet));
  auto new_rf2f_h = HostRead<LO>(new_rf2f);
  for (LO i = 0; i < nnew_tet; ++i) {
    LO const is_invalid = tet_invalid_h[i];
    if (is_invalid > 0) {
      //printf("writing file for tet %d\n", i);
      for (LO j=0; j<4; ++j) {
        LO const f = new_rf2f_h[i*4 + j];
        face_crvVis[f] = 1;
      }
      new_mesh->set_tag<LO>(2, "face_crvVis", Read<LO>(face_crvVis.write()));

      auto mesh_invalids = Mesh(new_mesh->comm()->library());
      mesh_invalids.set_comm(new_mesh->comm());
      build_cubic_cavities_3d(new_mesh, &mesh_invalids, 50);//curveVtk
      std::string vtuPath = "../omega_h/meshes/invalid_tet_";
      vtuPath += std::to_string(i);
      vtuPath += ".vtu";
      vtk::write_simplex_connectivity(vtuPath.c_str(), &mesh_invalids, 2);

      for (LO j=0; j<4; ++j) {
        LO const f = new_rf2f_h[i*4 + j];
        face_crvVis[f] = -1;
      }
      new_mesh->set_tag<LO>(2, "face_crvVis", Read<LO>(face_crvVis.write()));
    }
  }
  new_mesh->remove_tag(2, "face_crvVis");

  return;
}

void correct_curved_edges(Mesh *new_mesh) {

  auto const edge_crv2bdry_dim = new_mesh->get_array<I8>(1, "edge_crv2bdry_dim");
  auto const newedge_gid = new_mesh->get_array<LO>(1, "class_id");
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

  Write<LO> face_vis (nnew_face, -1);
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
          printf(" invalid tet %d for edge %d curved to g_dim %d, g_id %d with code %d\n"
              , adj_tet, i, edge_crv2bdry_dim[i], newedge_gid[i], has_invalid_tet);
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
  //new_mesh->set_tag<LO>(2, "face_crvVis", Read<LO>(face_vis));
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
