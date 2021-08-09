#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

Real intpow(const Real b, const LO e) {
  switch (e) {
  case 0: return 1.0;
  case 1: return b;
  case 2: return b*b;
  case 3: return b*b*b;
  case 4: return b*b*b*b;
  case 5: return b*b*b*b*b;
  case 6: return b*b*b*b*b*b;
  default:
    return intpow(b, e-6) * intpow(b, 6);
  }
}

Real B0(Real u) {
  return intpow(1-u, 3);
}

Real B1(Real u) {
  return 3*u*intpow(1-u, 2);
}

Real B2(Real u) {
  return 3*(1-u)*intpow(u, 2);
}

Real B3(Real u) {
  return intpow(u, 3);
}

void elevate_curve_order_2to3(Mesh* mesh) {

  I8 new_order = 3;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->set_max_order(new_order);
  auto n_new_pts = mesh->n_internal_ctrlPts(1);
  Write<Real> new_pts(nedge*n_new_pts*dim, 0.0);
  Write<Real> c1(dim, 0.0);
  Write<Real> c2(dim, 0.0);
  auto calc_pts = OMEGA_H_LAMBDA (LO i) {
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
  parallel_for(nedge, calc_pts);

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

void elevate_curve_order_3to4(Mesh* mesh) {

  I8 new_order = 4;
  auto old_ctrl_pts = mesh->get_ctrlPts(1);
  auto old_n_ctrl_pts = mesh->n_internal_ctrlPts(1);
  auto coords = mesh->coords();
  auto nedge = mesh->nedges();
  auto dim = mesh->dim();
  auto ev2v = mesh->get_adj(1, 0).ab2b;

  mesh->set_max_order(new_order);
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
              (3.0/4.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d];
      c2[d] = (2.0/4.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + d] +
              (2.0/4.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d];
      c3[d] = (3.0/4.0)*old_ctrl_pts[i*old_n_ctrl_pts*dim + 1*dim + d] +
              (1.0/4.0)*coords[v1*dim + d];
      new_pts[i*n_new_pts*dim + d] = c1[d];
      new_pts[i*n_new_pts*dim + 1*dim + d] = c2[d];
      new_pts[i*n_new_pts*dim + 2*dim + d] = c3[d];
    }
  };
  parallel_for(nedge, calc_pts);

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
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

  mesh->set_max_order(new_order);
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
  parallel_for(nedge, calc_pts);

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

  mesh->set_max_order(new_order);
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
  parallel_for(nedge, calc_pts);

  mesh->set_tag_for_ctrlPts(1, Reals(new_pts));
  return;
}

LO binomial(LO n, LO i) {

  i = std::min(n-i,i);

  if (i == 0)
    return 1;
  if (i == 1)
    return n;

  static int const bn4[1] = {6};
  static int const bn5[1] = {10};
  static int const bn6[2] = {15,20};
  static int const bn7[2] = {21,35};
  static int const bn8[3] = {28,56,70};
  static int const bn9[3] = {36,84,126};
  static int const bn10[4] = {45,120,210,252};
  static int const bn11[4] = {55,165,330,462};
  static int const bn12[5] = {66,220,495,792,924};
  static int const bn13[5] = {78,286,715,1287,1716};
  static int const bn14[6] = {91,364,1001,2002,3003,3432};
  static int const bn15[6] = {105,455,1365,3003,5005,6435};
  static int const bn16[7] = {120,560,1820,4368,8008,11440,12870};
  static int const bn17[7] = {136,680,2380,6188,12376,19448,24310};
  static int const bn18[8] = {153,816,3060,8568,18564,31824,43758,48620};
  static int const bn19[8] = {171,969,3876,11628,27132,50388,75582,92378};
  static int const bn20[9] = {190,1140,4845,15504,38760,77520,125970,167960,
      184756};
  static int const bn21[9] = {210,1330,5985,20349,54264,116280,203490,293930,
      352716};
  static int const bn22[10] = {231,1540,7315,26334,74613,170544,319770,497420,
      646646,705432};
  static int const bn23[10] = {253,1771,8855,33649,100947,245157,490314,817190,
      1144066,1352078};
  static int const bn24[11] = {276,2024,10626,42504,134596,346104,735471,1307504,
      1961256,2496144,2704156};
  static int const bn25[11] = {300,2300,12650,53130,177100,480700,1081575,2042975,
      3268760,4457400,5200300};
  static int const bn26[12] = {325,2600,14950,65780,230230,657800,1562275,3124550,
      5311735,7726160,9657700,10400600};
  static int const bn27[12] = {351,2925,17550,80730,296010,888030,2220075,4686825,
      8436285,13037895,17383860,20058300};
  static int const bn28[13] = {378,3276,20475,98280,376740,1184040,3108105,6906900,
      13123110,21474180,30421755,37442160,40116600};
  static int const bn29[13] = {406,3654,23751,118755,475020,1560780,4292145,10015005,
      20030010,34597290,51895935,67863915,77558760};
  static int const bn30[14] = {435,4060,27405,142506,593775,2035800,5852925,14307150,
      30045015,54627300,86493225,119759850,145422675,155117520};
  static int const bn31[14] = {465,4495,31465,169911,736281,2629575,7888725,20160075,
      44352165,84672315,141120525,206253075,265182525,300540195};
  static int const bn32[15] = {496,4960,35960,201376,906192,3365856,10518300,28048800,
      64512240,129024480,225792840,347373600,471435600,565722720,601080390};
  static int const bn33[15] = {528,5456,40920,237336,1107568,4272048,13884156,38567100,
      92561040,193536720,354817320,573166440,818809200,1037158320,1166803110};

  static int const* const bnTable[34] = {0,0,0,0,bn4,bn5,bn6,bn7,bn8,
      bn9,bn10,bn11,bn12,bn13,bn14,bn15,bn16,bn17,bn18,bn19,bn20,bn21,bn22,bn23,
      bn24,bn25,bn26,bn27,bn28,bn29,bn30,bn31,bn32,bn33};

  return bnTable[n][i-2];

}

OMEGA_H_DEVICE LO trinomial(int n, int i, int j) {
  return binomial (n, i) * binomial(n-i, j);
}

OMEGA_H_DEVICE LO quadnomial(int n, int i, int j, int k) {
  return binomial(n, i) * binomial(n-i, j) * binomial(n-i-j, k);
}

OMEGA_H_DEVICE Real Bij(const int i, const int j, const double u,
                        const double v) {
  return intpow(u,i) * intpow(v,j);
}

OMEGA_H_DEVICE Real Bij(const int ij[], const double xi[]) {
  return Bij(ij[0], ij[1], xi[0], xi[1]);
}

OMEGA_H_DEVICE Real Bijk(const int ijk[], const double xi[]) {
  return Bijk(ijk[0], ijk[1], ijk[2], xi[0], xi[1], xi[2]);
}

OMEGA_H_DEVICE Real Bijkl(const int ijkl[], const double xi[]) {
  return Bijkl(ijkl[0], ijkl[1], ijkl[2], ijkl[3], xi[0], xi[1], xi[2], xi[3]);
}

static unsigned const b2_0_0[1] = {2};
static unsigned const* const b2_0[1] = {b2_0_0};
static unsigned const b2_1_0[2] = {2,1};
static unsigned const b2_1_1[1] = {0};
static unsigned const* const b2_1[2] = {b2_1_0,b2_1_1};
static unsigned const b2_2_0[3] = {2,4,1};
static unsigned const b2_2_1[2] = {5,3};
static unsigned const b2_2_2[1] = {0};
static unsigned const* const b2_2[3] = {b2_2_0,b2_2_1,b2_2_2};
static unsigned const b2_3_0[4] = {2,6,5,1};
static unsigned const b2_3_1[3] = {7,9,4};
static unsigned const b2_3_2[2] = {8,3};
static unsigned const b2_3_3[1] = {0};
static unsigned const* const b2_3[4] = {b2_3_0,b2_3_1,b2_3_2,b2_3_3};
static unsigned const b2_4_0[5] = {2,8,7,6,1};
static unsigned const b2_4_1[4] = {9,14,13,5};
static unsigned const b2_4_2[3] = {10,12,4};
static unsigned const b2_4_3[2] = {11,3};
static unsigned const b2_4_4[1] = {0};
static unsigned const* const b2_4[5] = {b2_4_0,b2_4_1,b2_4_2,b2_4_3,b2_4_4};
static unsigned const b2_5_0[6] = {2,10,9,8,7,1};
static unsigned const b2_5_1[5] = {11,20,19,17,6};
static unsigned const b2_5_2[4] = {12,18,16,5};
static unsigned const b2_5_3[3] = {13,15,4};
static unsigned const b2_5_4[2] = {14,3};
static unsigned const b2_5_5[1] = {0};
static unsigned const* const b2_5[6] =
{b2_5_0,b2_5_1,b2_5_2,b2_5_3,b2_5_4,b2_5_5};
static unsigned const b2_6_0[7] = {2,12,11,10,9,8,1};
static unsigned const b2_6_1[6] = {13,27,26,24,21,7};
static unsigned const b2_6_2[5] = {14,25,23,20,6};
static unsigned const b2_6_3[4] = {15,22,19,5};
static unsigned const b2_6_4[3] = {16,18,4};
static unsigned const b2_6_5[2] = {17,3};
static unsigned const b2_6_6[1] = {0};
static unsigned const* const b2_6[7] =
{b2_6_0,b2_6_1,b2_6_2,b2_6_3,b2_6_4,b2_6_5,b2_6_6};
static unsigned const b2_7_0[8] = {2,14,13,12,11,10,9,1};
static unsigned const b2_7_1[7] = {15,35,34,32,29,25,8};
static unsigned const b2_7_2[6] = {16,33,31,28,24,7};
static unsigned const b2_7_3[5] = {17,30,27,23,6};
static unsigned const b2_7_4[4] = {18,26,22,5};
static unsigned const b2_7_5[3] = {19,21,4};
static unsigned const b2_7_6[2] = {20,3};
static unsigned const b2_7_7[1] = {0};
static unsigned const* const b2_7[8] =
{b2_7_0,b2_7_1,b2_7_2,b2_7_3,b2_7_4,b2_7_5,b2_7_6,b2_7_7};
static unsigned const b2_8_0[9] = {2,16,15,14,13,12,11,10,1};
static unsigned const b2_8_1[8] = {17,44,43,41,38,34,29,9};
static unsigned const b2_8_2[7] = {18,42,40,37,33,28,8};
static unsigned const b2_8_3[6] = {19,39,36,32,27,7};
static unsigned const b2_8_4[5] = {20,35,31,26,6};
static unsigned const b2_8_5[4] = {21,30,25,5};
static unsigned const b2_8_6[3] = {22,24,4};
static unsigned const b2_8_7[2] = {23,3};
static unsigned const b2_8_8[1] = {0};
static unsigned const* const b2_8[9] =
{b2_8_0,b2_8_1,b2_8_2,b2_8_3,b2_8_4,b2_8_5,b2_8_6,b2_8_7,b2_8_8};
static unsigned const b2_9_0[10] = {2,18,17,16,15,14,13,12,11,1};
static unsigned const b2_9_1[9] = {19,54,53,51,48,44,39,33,10};
static unsigned const b2_9_2[8] = {20,52,50,47,43,38,32,9};
static unsigned const b2_9_3[7] = {21,49,46,42,37,31,8};
static unsigned const b2_9_4[6] = {22,45,41,36,30,7};
static unsigned const b2_9_5[5] = {23,40,35,29,6};
static unsigned const b2_9_6[4] = {24,34,28,5};
static unsigned const b2_9_7[3] = {25,27,4};
static unsigned const b2_9_8[2] = {26,3};
static unsigned const b2_9_9[1] = {0};
static unsigned const* const b2_9[10] =
{b2_9_0,b2_9_1,b2_9_2,b2_9_3,b2_9_4,b2_9_5,b2_9_6,b2_9_7,b2_9_8,b2_9_9};
static unsigned const b2_10_0[11] = {2,20,19,18,17,16,15,14,13,12,1};
static unsigned const b2_10_1[10] = {21,65,64,62,59,55,50,44,37,11};
static unsigned const b2_10_2[9] = {22,63,61,58,54,49,43,36,10};
static unsigned const b2_10_3[8] = {23,60,57,53,48,42,35,9};
static unsigned const b2_10_4[7] = {24,56,52,47,41,34,8};
static unsigned const b2_10_5[6] = {25,51,46,40,33,7};
static unsigned const b2_10_6[5] = {26,45,39,32,6};
static unsigned const b2_10_7[4] = {27,38,31,5};
static unsigned const b2_10_8[3] = {28,30,4};
static unsigned const b2_10_9[2] = {29,3};
static unsigned const b2_10_10[1] = {0};
static unsigned const* const b2_10[11] =
{b2_10_0,b2_10_1,b2_10_2,b2_10_3,b2_10_4,b2_10_5,b2_10_6,b2_10_7,b2_10_8,
    b2_10_9,b2_10_10};

unsigned const* const* const b2[11] =
{b2_0,b2_1,b2_2,b2_3,b2_4,b2_5,b2_6,b2_7,b2_8,b2_9,b2_10};

static unsigned const b3_0_00[1] = {3};
static unsigned const* const b3_0_0[1] = {b3_0_00};
static unsigned const* const* const b3_0[1] = {b3_0_0};

static unsigned const b3_1_00[2] = {3,2};
static unsigned const b3_1_01[1] = {1};
static unsigned const b3_1_10[1] = {0};
static unsigned const* const b3_1_0[2] = {b3_1_00,b3_1_01};
static unsigned const* const b3_1_1[1] = {b3_1_10};
static unsigned const* const* const b3_1[2] = {b3_1_0,b3_1_1};

static unsigned const b3_2_00[3] = {3,9,2};
static unsigned const b3_2_01[2] = {8,5};
static unsigned const b3_2_02[1] = {1};
static unsigned const b3_2_10[2] = {7,6};
static unsigned const b3_2_11[1] = {4};
static unsigned const b3_2_20[1] = {0};
static unsigned const* const b3_2_0[3] = {b3_2_00,b3_2_01,b3_2_02};
static unsigned const* const b3_2_1[2] = {b3_2_10,b3_2_11};
static unsigned const* const b3_2_2[1] = {b3_2_20};
static unsigned const* const* const b3_2[3] = {b3_2_0,b3_2_1,b3_2_2};

static unsigned const b3_3_00[4] = {3,15,14,2};
static unsigned const b3_3_01[3] = {13,18,7};
static unsigned const b3_3_02[2] = {12,6};
static unsigned const b3_3_03[1] = {1};
static unsigned const b3_3_10[3] = {11,19,8};
static unsigned const b3_3_11[2] = {17,16};
static unsigned const b3_3_12[1] = {5};
static unsigned const b3_3_20[2] = {10,9};
static unsigned const b3_3_21[1] = {4};
static unsigned const b3_3_30[1] = {0};
static unsigned const* const b3_3_0[4] = {b3_3_00,b3_3_01,b3_3_02,b3_3_03};
static unsigned const* const b3_3_1[3] = {b3_3_10,b3_3_11,b3_3_12};
static unsigned const* const b3_3_2[2] = {b3_3_20,b3_3_21};
static unsigned const* const b3_3_3[1] = {b3_3_30};
static unsigned const* const* const b3_3[4] =
  {b3_3_0,b3_3_1,b3_3_2,b3_3_3};

static unsigned const b3_4_00[5] = {3,21,20,19,2};
static unsigned const b3_4_01[4] = {18,30,29,9};
static unsigned const b3_4_02[3] = {17,28,8};
static unsigned const b3_4_03[2] = {16,7};
static unsigned const b3_4_04[1] = {1};
static unsigned const b3_4_10[4] = {15,33,32,10};
static unsigned const b3_4_11[3] = {27,34,24};
static unsigned const b3_4_12[2] = {26,23};
static unsigned const b3_4_13[1] = {6};
static unsigned const b3_4_20[3] = {14,31,11};
static unsigned const b3_4_21[2] = {25,22};
static unsigned const b3_4_22[1] = {5};
static unsigned const b3_4_30[2] = {13,12};
static unsigned const b3_4_31[1] = {4};
static unsigned const b3_4_40[1] = {0};
static unsigned const* const b3_4_0[5] =
  {b3_4_00,b3_4_01,b3_4_02,b3_4_03,b3_4_04};
static unsigned const* const b3_4_1[4] = {b3_4_10,b3_4_11,b3_4_12,b3_4_13};
static unsigned const* const b3_4_2[3] = {b3_4_20,b3_4_21,b3_4_22};
static unsigned const* const b3_4_3[2] = {b3_4_30,b3_4_31};
static unsigned const* const b3_4_4[1] = {b3_4_40};
static unsigned const* const* const b3_4[5] =
  {b3_4_0,b3_4_1,b3_4_2,b3_4_3,b3_4_4};

unsigned const* const* const* const b3[5] =
{b3_0,b3_1,b3_2,b3_3,b3_4};

OMEGA_H_DEVICE LO computeTriNodeIndex(int P, int i, int j) {
  int k = P-i-j;
  if(i == P) return 0;
  if(j == P) return 1;
  if(k == P) return 2;
  if(k == 0) return 2+j; // 0-1
  if(i == 0) return 2+(P-1)+k; // 1-2
  if(j == 0) return 2+(P-1)*2+i; // 2-0
  return k*(P-1)-k*(k-1)/2+j+2*P;
}

OMEGA_H_DEVICE LO computeTetNodeIndex(int P, int i, int j, int k) {
  int l = P-i-j-k;
  if(i == P) return 0;
  if(j == P) return 1;
  if(k == P) return 2;
  if(l == P) return 3;
  if(k == 0 && l == 0) return 3+j; // 0-1
  if(i == 0 && l == 0) return 3+(P-1)+k; // 1-2
  if(j == 0 && l == 0) return 3+2*(P-1)+i; // 2-0
  if(j == 0 && k == 0) return 3+3*(P-1)+l; // 0-3
  if(i == 0 && k == 0) return 3+4*(P-1)+l; // 1-3
  if(i == 0 && j == 0) return 3+5*(P-1)+l;// 2-3
  if(l == 0) return k*(P-1)-k*(k-1)/2+j+5*P-2; // 0-1-2
  if(k == 0) return l*(P-1)-l*(l-1)/2+j+5*P-2+(P-2)*(P-1)/2; // 0-1-3
  if(i == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1);// 1-2-3
  if(j == 0) return l*(P-1)-l*(l-1)/2+k+5*P-2+(P-2)*(P-1)*3/2; // 0-2-3
  return i-P-((i-P+1)*(i-P+2)*(i-P+3))/6+l*(P-1-i)-l*(l-1)/2+k+2*P*P+2;
}

OMEGA_H_DEVICE LO getTriNodeIndex(int P, int i, int j) {
  // use a table if its small, otherwise dynamically generate it on the fly
  if(P <= 10)
    return b2[P][i][j];
  else
    return computeTriNodeIndex(P,i,j);
}

OMEGA_H_DEVICE LO getTetNodeIndex(int P, int i, int j, int k) {
  if(P <= 4)
    return b3[P][i][j][k];
  else
    return computeTetNodeIndex(P,i,j,k);
}

/*
OMEGA_H_DEVICE static void bezierCurve(I8 P, Reals xi, Write<Real> values) {
  double t = 0.5 * (xi[0] + 1.0);
  for (I8 i = 1; i < P; ++i) {
    values[i+1] = binomial(P,i) * Bij(P-i, i, 1.0-t, t);
  }
  values[0] = intpow(1-t, P);
  values[1] = intpow(t, P);
}

OMEGA_H_DEVICE static void bezierTriangle(I8 P, Reals xi, Write<Real> values) {
  double xii[3] = {1.0-xi[0]-xi[1], xi[0], xi[1]};
  for (int i = 0; i < P+1; ++i) {
    for (int j = 0; j < P+1-i; ++j) {
      values[getTriNodeIndex(P, i, j)] =
        trinomial(P, i, j)*Bijk(i, j, P-i-j, xii[0], xii[1], xii[2]);
    }
  }
}

OMEGA_H_DEVICE static void bezierTet(I8 P, Reals xi, Write<Real> values) {
  double xii[4] = {1.0-xi[0]-xi[1]-xi[2], xi[0], xi[1], xi[2]};
  for (I8 i = 0; i < 4; ++i) {
    values[i] = intpow(xii[i], P);
  }

  int nE = P - 1;

  for (int a = 0; a < 6; ++a) {
    for (int b = 0; b < nE; ++b) { // edge nodes
      values[4+a*nE+b] = binomial(P, b+1)
        *Bij(P-b-1, b+1, xii[simplex_down_template(3, 1, a, 0)],
             xii[simplex_down_template(3, 1, a, 1)]);
    }
  }

  // face 0, l = 0
  for (int i = 1; i <= P-1; ++i) {
    for (int j = 1; j <= P-1-i; ++j) {
      values[computeTetNodeIndex(P, i, j, P-i-j)] = trinomial(P, i, j)
        *Bijk(i, j, P-i-j, xii[0], xii[1], xii[2]);
    }
  }
  // face 1, k = 0
  for (int i = 1; i <= P-1; ++i) {
    for (int j = 1; j <= P-1-i; ++j) {
      values[computeTetNodeIndex(P, i, j, 0)] = trinomial(P, i, j)
        *Bijk(i, j, P-i-j, xii[0], xii[1], xii[3]);
    }
  }
  // face 2, i = 0
  for (int j = 1; j <= P-1; ++j) {
    for (int k = 1; k <= P-1-j; ++k) {
      values[computeTetNodeIndex(P, 0, j, k)] = trinomial(P, j, k)
        *Bijk(j, k, P-j-k, xii[1], xii[2], xii[3]);
    }
  }
  // face 3, j = 0
  for (int i = 1; i <= P-1; ++i) {
    for (int k = 1; k <= P-1-i; ++k) {
      values[computeTetNodeIndex(P, i, 0, k)] = trinomial(P, i, k)
        *Bijk(i, k, P-i-k, xii[0], xii[2], xii[3]);
    }
  }

  // internal nodes
  for (int i = 1; i <= P-1; ++i) {
    for (int j = 1; j <= P-1-i; ++j) {
      for (int k = 1; k <= P-1-i-j; ++k) {
        values[computeTetNodeIndex(P, i, j, k)] = quadnomial(P, i, j, k)*
          Bijkl(i, j, k, P-i-j-k, xii[0], xii[1], xii[2], xii[3]);
      }
    }
  }
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
