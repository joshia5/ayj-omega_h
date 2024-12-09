#include <Omega_h_build.hpp>
#include <Omega_h_swap.hpp>
#include <Omega_h_library.hpp>
//#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_beziers.hpp>
#include <Omega_h_bezier_interp.hpp>
//#include <Omega_h_for.hpp>
#include <Omega_h_curve_coarsen.hpp>
#include <Omega_h_curve_validity_3d.hpp>
#include <LBFGS.hpp>

using namespace Omega_h;

#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>

using Eigen::VectorXf;
using Eigen::MatrixXf;
using namespace LBFGSpp;

class Rosenbrock
{
  private:
    int n;
  public:
    Rosenbrock(int n_) : n(n_) {}
    float operator()(const VectorXf& x, VectorXf& grad)
    {   
      float fx = 0.0;
      for(int i = 0; i < n; i += 2)
      {
        float t1 = 1.0 - x[i];
        float t2 = 10 * (x[i + 1] - x[i] * x[i]);
        grad[i + 1] = 20 * t2; 
        grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
        fx += t1 * t1 + t2 * t2; 
      }
      return fx; 
    }   
};

void test_rosenbrock() {

  const int n = 10; 
  LBFGSParam<float> param;
  LBFGSSolver<float> solver(param);
  Rosenbrock fun(n);

  VectorXf x = VectorXf::Zero(n);
  float fx; 
  int niter = solver.minimize(fun, x, fx);

  std::cout << niter << " iterations" << std::endl;
  std::cout << "x = \n" << x.transpose() << std::endl;
  std::cout << "f(x) = " << fx << std::endl;
  std::cout << "grad = " << solver.final_grad().transpose() << std::endl;
  std::cout << "||grad|| = " << solver.final_grad_norm() << std::endl;

  return;
}

void test_annulus_optim(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-120d-4.sms",
                            "/lore/joshia5/Models/curved/annulus-120cut.smd", comm);
 
  //auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus-8.sms",
    //                        "/lore/joshia5/Models/curved/annulus-8.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());
  auto const dim = mesh.dim();

  auto cubic_2tri = Mesh(lib);
  cubic_2tri.set_comm(comm);
  cubic_2tri.set_parting(OMEGA_H_ELEM_BASED);
  cubic_2tri.set_dim(dim);
  cubic_2tri.set_family(OMEGA_H_SIMPLEX);
  LO const numVtx=4;
  LO const numEdge=5;
  LO const numFace=2;
  cubic_2tri.set_verts(numVtx);

  HostWrite<Real> host_coords(numVtx*dim);
  host_coords = {
    0,0.25,
    0,0.5,
    -0.433013,-0.25,
    -0.216506,-0.125
                };
  HostWrite<LO> class_id_0(numVtx);
  class_id_0 = {12,11,16,14};
  HostWrite<I8> class_dim_0(numVtx);// all 0s
  class_dim_0 = {0,0,0,0};
  cubic_2tri.add_coords(Reals(host_coords.write()));
  cubic_2tri.add_tag<ClassId>(0, "class_id", 1,
                           Read<LO>(class_id_0.write()));
  cubic_2tri.add_tag<I8>(0, "class_dim", 1,
                      Read<I8>(class_dim_0.write()));

  HostWrite<LO> ev2v(numEdge*2);
  ev2v = {0,1,  1,2,  2,3,  3,0, 1,3};
  cubic_2tri.set_ents(1, Adj(Read<LO>(ev2v.write())));
  HostWrite<LO> class_id_1(numEdge);
  class_id_1 = {13,7,17,15,2};
  HostWrite<I8> class_dim_1(numEdge);// all 0s
  class_dim_1 = {1,1,1,1,2};
  cubic_2tri.add_tag<ClassId>(1, "class_id", 1,
                           Read<LO>(class_id_1.write()));
  cubic_2tri.add_tag<I8>(1, "class_dim", 1,
                      Read<I8>(class_dim_1.write()));
  
  HostWrite<LO> fv2v(numFace*3);
  fv2v = {0,1,3,  2,3,1};
  //cubic_2tri.set_down(2, 0, LOs(fv2v.write()));//cant set if using set_ents
  Adj edge2vert = cubic_2tri.get_adj(1, 0);
  Adj vert2edge = cubic_2tri.ask_up(0, 1);
  Adj f2e = reflect_down(LOs(fv2v.write()), edge2vert.ab2b, vert2edge,
      OMEGA_H_SIMPLEX, 2, 1);
  cubic_2tri.set_ents(2, f2e);
  printf("ok5\n");
  HostWrite<LO> class_id_2(numFace);
  class_id_2 = {2,2};
  HostWrite<I8> class_dim_2(numFace);// all 0s
  class_dim_2 = {2,2};
  cubic_2tri.add_tag<ClassId>(2, "class_id", 1,
                           Read<LO>(class_id_2.write()));
  cubic_2tri.add_tag<I8>(2, "class_dim", 1,
                      Read<I8>(class_dim_2.write()));
  for (LO i = 0; i <= cubic_2tri.dim(); ++i) {
    if (!cubic_2tri.has_tag(i, "global")) {
      cubic_2tri.add_tag(i, "global", 1, Omega_h::GOs(cubic_2tri.nents(i), 0, 1));
    }
  }
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-120-2tri_full.vtk", &cubic_2tri);
  writer.write();
  cubic_2tri.set_curved(1);
  cubic_2tri.set_max_order(3);
  cubic_2tri.add_tags_for_ctrlPts();
  cubic_2tri.add_tag<Real>(0, "bezier_pts", cubic_2tri.dim(), cubic_2tri.coords());
  printf("making 2tri curve mesh\n");

  HostWrite<Real> edgePt_coords(numEdge*dim*2);
  // init straight sided
  for (LO i=0; i<numEdge; ++i) {
    LO const v0 = ev2v[i*2 + 0];
    LO const v1 = ev2v[i*2 + 1];
    Real len = 0.;
    for (LO d=0; d<dim; ++d) {
      len += std::pow((host_coords[v1*dim + d] - host_coords[v0*dim + d]),2);
    }
    len = std::pow(len, 0.5);
    for (LO d=0; d<dim; ++d) {
      edgePt_coords[i*dim*2 + d] = host_coords[v0*dim + d] + 
        (host_coords[v1*dim + d] - host_coords[v0*dim + d])/3.;
      edgePt_coords[i*dim*2 + dim + d] = host_coords[v0*dim + d] + 
        (host_coords[v1*dim + d] - host_coords[v0*dim + d])*2./3.;
    }
    Vector<2> c0,c3,p1,p2;
    c0[0] = host_coords[v0*dim + 0]; 
    c0[1] = host_coords[v0*dim + 1]; 
    c3[0] = host_coords[v1*dim + 0]; 
    c3[1] = host_coords[v1*dim + 1]; 
    if (i == 1) {
      /*
      p1[0] = 0.5*std::cos(PI*130./180.);
      p1[1] = 0.5*std::sin(PI*130./180.);
      p2[0] = 0.5*std::cos(PI*170./180.);
      p2[1] = 0.5*std::sin(PI*170./180.);

      auto const c1_c2 = curve_interpToCtrl_pts_2d(3, c0, c3, p1, p2);
      edgePt_coords[i*dim*2 + 0] = c1_c2[0];
      edgePt_coords[i*dim*2 + 1] = c1_c2[1];
      edgePt_coords[i*dim*2 + dim + 0] = c1_c2[2];
      edgePt_coords[i*dim*2 + dim + 1] = c1_c2[3];
      */

      edgePt_coords[i*dim*2 + 0] = 0.5*std::cos(PI*130./180.);
      edgePt_coords[i*dim*2 + 1] = 0.5*std::sin(PI*130./180.);
      
      edgePt_coords[i*dim*2 + dim + 0] = 0.5*std::cos(PI*170./180.);
      edgePt_coords[i*dim*2 + dim + 1] = 0.5*std::sin(PI*170./180.);
    }
    if (i == 3) {
      edgePt_coords[i*dim*2 + 0] = 0.25*std::cos(PI*170./180.);
      edgePt_coords[i*dim*2 + 1] = 0.25*std::sin(PI*170./180.);
      
      edgePt_coords[i*dim*2 + dim + 0] = 0.25*std::cos(PI*130./180.);
      edgePt_coords[i*dim*2 + dim + 1] = 0.25*std::sin(PI*130./180.);
    }
  }
  HostWrite<Real> facePt_coords(numFace*dim);
  for (LO i=0; i<numFace; ++i) {
    LO const v0 = fv2v[i*2 + 0];
    LO const v1 = fv2v[i*2 + 1];
    LO const v2 = fv2v[i*2 + 2];
    for (LO d=0; d<dim; ++d) {
      facePt_coords[i*dim + d] = (host_coords[v0*dim + d] + 
        host_coords[v1*dim + d] + host_coords[v2*dim + d])/3.;
    }
  }
  cubic_2tri.set_tag_for_ctrlPts(1, Reals(edgePt_coords.write()));
  cubic_2tri.set_tag_for_ctrlPts(2, Reals(facePt_coords.write()));
  
  auto opts = AdaptOpts(&cubic_2tri);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  cubic_2tri.add_tag<Real>(VERT, "metric", 1);
  cubic_2tri.set_tag(VERT, "metric", Reals(cubic_2tri.nverts(), 1));
  auto valid_tris_bef = checkValidity_2d(&cubic_2tri, LOs(cubic_2tri.nfaces(), 0, 1), 2);

  /*
  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
  auto valid_tris_bef = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  */
  askWorstQuality_2d(&cubic_2tri, LOs(cubic_2tri.nfaces(), 0, 1), 2);
  //auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&cubic_2tri, opts);
    LBFGS *l;
    ObjFunction *objF;
    std::vector<double> x0;
    for (LO i = 0; i < 4*2; ++i) {
      x0.push_back(0.);
    }
    l = new LBFGS (0.001, 5, x0, objF);
    auto lval = l->run(&cubic_2tri);
  }
/*
  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
*/
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  //form vandermonde matrix using bernsteins
  //form ctrl pt vector
  //for interp pt vector

  return;
}

void test_annulus3d_swap(Library *lib) {
  auto comm = lib->world();

  auto mesh = meshsim::read("/lore/joshia5/Meshes/curved/annulus3d-24.sms",
                            "/lore/joshia5/Models/curved/annulus3d.smd", comm);
 
  calc_quad_ctrlPts_from_interpPts(&mesh);
  elevate_curve_order_2to3(&mesh);
  mesh.add_tag<Real>(0, "bezier_pts", mesh.dim(), mesh.coords());

  auto wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_3d(&mesh, &wireframe_mesh, 4);
  std::string vtuPath =
    "/lore/joshia5/Meshes/curved/annulus3d_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  auto cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_3d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus3d.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);

  auto opts = AdaptOpts(&mesh);
  opts.should_coarsen = false;
  opts.should_coarsen_slivers = false;
  opts.should_refine = false;
  opts.should_filter_invalids = false;
  opts.verbosity = EXTRA_STATS;
  opts.min_quality_desired = 0.99;
  opts.min_quality_allowed = 0.98;
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1));
  auto valid_tris_bef = checkValidity_3d(&mesh, LOs(mesh.nregions(), 0, 1));
  auto qual = calc_crvQuality_3d(&mesh);
  /*
  //auto quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  for (LO adapt_itr = 0; adapt_itr < 1; ++adapt_itr) {
    fprintf(stderr, "itr %d\n", adapt_itr);
    swap_edges(&mesh, opts);
  }

  wireframe_mesh = Mesh(lib);
  wireframe_mesh.set_comm(comm);
  build_cubic_wireframe_2d(&mesh, &wireframe_mesh, 4);
  vtuPath =
    "/lore/joshia5/Meshes/curved/annulus-8-swap_wire.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &wireframe_mesh, 1);
  cubic_curveVtk_mesh = Mesh(lib);
  cubic_curveVtk_mesh.set_comm(comm);
  build_cubic_curveVtk_2d(&mesh, &cubic_curveVtk_mesh, 4);
  vtuPath = "/lore/joshia5/Meshes/curved/annulus-8-swap.vtu";
  vtk::write_simplex_connectivity(vtuPath.c_str(), &cubic_curveVtk_mesh, 2);
  auto valid_tris_aft = checkValidity_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);
  //quals = askQuality_2d(&mesh, LOs(mesh.nfaces(), 0, 1), 2);

  fprintf(stderr, "finish swapping annulus-3d\n");
  vtk::FullWriter writer;
  writer = vtk::FullWriter(
      "/lore/joshia5/Meshes/curved/annulus-8-swap_full.vtk", &mesh);
  writer.write();
  */
  //mesh.ask_qualities();

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);

  //test_annulus_optim(&lib); //2d

  //test_annulus3d_swap(&lib); //2d
  //
  
  test_rosenbrock();
  return 0;
}
