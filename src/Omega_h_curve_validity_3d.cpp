#include "Omega_h_curve_validity_3d.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

LOs checkValidity_3d(Mesh *mesh, LOs new_tets) {
  auto rv2v = mesh->ask_down(3, 0).ab2b;
  auto re2e = mesh->ask_down(3, 1).ab2b;
  auto rf2f = mesh->get_adj(3, 2).ab2b;
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);
  auto faceCtrlPts = mesh->get_ctrlPts(2);
  auto order = mesh->get_max_order();
  OMEGA_H_CHECK(order <= 3);

  Write<LO> is_invalid(new_tets.size(), -1);

  auto check_validity = OMEGA_H_LAMBDA (LO n) {
    auto tet = new_tets[n];

    Few<Real, 60> tet_pts = collect_tet_pts(order, tet, ev2v, rv2v, vertCtrlPts
        , edgeCtrlPts, faceCtrlPts, re2e, rf2f);
    
    Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

    is_invalid[n] = checkMinJacDet_3d(nodes_det, order);
    if (is_invalid[n] > 0) printf("invalid tet %d\n", tet);
  };
  parallel_for(new_tets.size(), std::move(check_validity));

  return LOs(is_invalid);
}

}
