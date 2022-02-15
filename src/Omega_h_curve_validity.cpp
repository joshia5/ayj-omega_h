#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_few.hpp"
#include "Omega_h_scalar.hpp"

namespace Omega_h {

LO binomial(int n, int i) {

  i = std::min(n-i,i);

  if (i == 0)
    return 1;
  if (i == 1)
    return n;

  static LO const bn4[1] = {6};
  static LO const bn5[1] = {10};
  static LO const bn6[2] = {15,20};
  static LO const bn7[2] = {21,35};
  static LO const bn8[3] = {28,56,70};
  static LO const bn9[3] = {36,84,126};
  static LO const bn10[4] = {45,120,210,252};
  static LO const bn11[4] = {55,165,330,462};
  static LO const bn12[5] = {66,220,495,792,924};
  static LO const bn13[5] = {78,286,715,1287,1716};
  static LO const bn14[6] = {91,364,1001,2002,3003,3432};
  static LO const bn15[6] = {105,455,1365,3003,5005,6435};
  static LO const bn16[7] = {120,560,1820,4368,8008,11440,12870};
  static LO const bn17[7] = {136,680,2380,6188,12376,19448,24310};
  static LO const bn18[8] = {153,816,3060,8568,18564,31824,43758,48620};
  static LO const bn19[8] = {171,969,3876,11628,27132,50388,75582,92378};
  static LO const bn20[9] = {190,1140,4845,15504,38760,77520,125970,167960,
    184756};
  static LO const bn21[9] = {210,1330,5985,20349,54264,116280,203490,293930,
    352716};
  static LO const bn22[10] = {231,1540,7315,26334,74613,170544,319770,497420,
    646646,705432};
  static LO const bn23[10] = {253,1771,8855,33649,100947,245157,490314,817190,
    1144066,1352078};
  static LO const bn24[11] = {276,2024,10626,42504,134596,346104,735471,1307504,
    1961256,2496144,2704156};
  static LO const bn25[11] = {300,2300,12650,53130,177100,480700,1081575,2042975,
    3268760,4457400,5200300};
  static LO const bn26[12] = {325,2600,14950,65780,230230,657800,1562275,3124550,
    5311735,7726160,9657700,10400600};
  static LO const bn27[12] = {351,2925,17550,80730,296010,888030,2220075,4686825,
    8436285,13037895,17383860,20058300};
  static LO const bn28[13] = {378,3276,20475,98280,376740,1184040,3108105,6906900,
    13123110,21474180,30421755,37442160,40116600};
  static LO const bn29[13] = {406,3654,23751,118755,475020,1560780,4292145,10015005,
    20030010,34597290,51895935,67863915,77558760};
  static LO const bn30[14] = {435,4060,27405,142506,593775,2035800,5852925,14307150,
    30045015,54627300,86493225,119759850,145422675,155117520};
  static LO const bn31[14] = {465,4495,31465,169911,736281,2629575,7888725,20160075,
    44352165,84672315,141120525,206253075,265182525,300540195};
  static LO const bn32[15] = {496,4960,35960,201376,906192,3365856,10518300,28048800,
    64512240,129024480,225792840,347373600,471435600,565722720,601080390};
  static LO const bn33[15] = {528,5456,40920,237336,1107568,4272048,13884156,38567100,
    92561040,193536720,354817320,573166440,818809200,1037158320,1166803110};

  static LO const* const bnTable[34] = {0,0,0,0,bn4,bn5,bn6,bn7,bn8,
    bn9,bn10,bn11,bn12,bn13,bn14,bn15,bn16,bn17,bn18,bn19,bn20,bn21,bn22,bn23,
    bn24,bn25,bn26,bn27,bn28,bn29,bn30,bn31,bn32,bn33};

  return bnTable[n][i-2];
}

LO trinomial(LO n, LO i, LO j) {
  return binomial(n,i)*binomial(n-i,j);
}

OMEGA_H_INLINE LO computeTriNodeIndex (LO P, LO i, LO j) {
  LO k = P-i-j;
  if (i == P) return 0;
  if (j == P) return 1;
  if (k == P) return 2;
  if (k == 0) return 2+j;
  if (i == 0) return 2+(P-1)+k;
  if (j == 0) return 2+(P-1)*2+i;
  return k*(P-1)-k*(k-1)/2+j+2*P;
}

OMEGA_H_INLINE LO getTriNodeIndex (LO P, LO i, LO j) {
  return computeTriNodeIndex(P,i,j);
}

template <Int n>
OMEGA_H_INLINE Real getTriPartialJacobianDet(Few<Real, 20> nodes,
  LO P, LO i1, LO j1, LO i2, LO j2) {
  LO p00 = getTriNodeIndex(P,i1+1,j1);
  LO p01 = getTriNodeIndex(P,i1,j1+1);
  LO p10 = getTriNodeIndex(P,i2+1,j2);
  LO p11 = getTriNodeIndex(P,i2,j2);
  return cross(get_vector<n>(nodes, p01) - get_vector<n>(nodes, p00),
               get_vector<n>(nodes, p11) - get_vector<n>(nodes, p10));
  /*the return will be following for 3d
  return cross(get_vector<n>(nodes, p01) - get_vector<n>(nodes, p00),
               get_vector<n>(nodes, p11) - get_vector<n>(nodes, p10))[2];
  */
}

OMEGA_H_INLINE Real Nijk(Few<Real, 20> nodes, LO d, LO I, LO J) {
  Real sum = 0.;
  LO CD = trinomial(2*(d-1), I, J);
  for (LO j1 = 0; j1 <= J; ++j1) {
    auto i1start = max2(0, I+J-j1-(d-1));
    auto i1end = min2(I, d-1-j1);
    for (LO i1 = i1start; i1 <= i1end; ++i1){
      sum += trinomial(d-1, i1, j1)*trinomial(d-1, I-i1, J-j1)
        *getTriPartialJacobianDet<2>(nodes, d, i1, j1, I-i1, J-j1);
    }
  }
  return sum*d*d/CD;
}

template <Int n>
OMEGA_H_INLINE Few<Real, n> getTriJacDetNodes(
    LO P, Few<Real, 20> const& elemNodes) {
  Few<Real, n> nodes;//n=15
  for (LO I = 0; I <= 2*(P-1); ++I) {
    for (LO J = 0; J <= 2*(P-1)-I; ++J) {
        OMEGA_H_CHECK(getTriNodeIndex(2*(P-1),I,J) < nodes.size());
        nodes[getTriNodeIndex(2*(P-1),I,J)] = Nijk(elemNodes,P,I,J);
    }
  }
  return nodes;
}

template<Int n>
OMEGA_H_INLINE LO checkMinJacDet(Few<Real, n> const& nodes) {
  // first 3 vertices
  Real minAcceptable = 0.0;
  for (LO i = 0; i < 3; ++i) {
    fprintf(stderr, "i %d Nijk[i] %f\n", i, nodes[i]);
    if (nodes[i] < minAcceptable) {
    //if (std::abs(nodes[i] - minAcceptable) > EPSILON) {
      return i+2;
    }
  }

  /*
  Real minJ = 0;
  for (LO edge = 0; edge < 3; ++edge) {
    for (LO i = 0; i < 2*(order-1)-1; ++i) {
      if (nodes[3+edge*(2*(order-1)-1)+i] < minAcceptable) {
        minJ = -1e10;
        // there is no point in doing much with edges if we dont have
        // elevation or subdivision or subdivision matrices
        if (minJ < minAcceptable){
          return 8+edge;
        }
      }
    }
  }

  for (LO i = 0; i < (2*order-3)*(2*order-4)/2; ++i) {
    if (nodes[6*(order-1)+i] < minAcceptable) {
      minJ = -1e10;
      if (minJ < minAcceptable) {
        return 14;
      }
    }
  }
  */
  return -1;
}

LOs checkValidity(Mesh *mesh, LOs new_tris, Int const mesh_dim) {
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

    auto nodes_det = getTriJacDetNodes<15>(order, tri_pts);

    is_invalid[n] = checkMinJacDet<15>(nodes_det);
    if (is_invalid[n] > 0) Omega_h_fail("invalid tri %d\n", tri);
  };
  parallel_for(new_tris.size(), std::move(check_validity));

  return LOs(is_invalid);
}

} //namespace
