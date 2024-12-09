#ifndef LBFGS_H
#define LBFGS_H

#include <iostream>
//#include <pcu_util.h>
#include <vector>

#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"

namespace Omega_h {

class ObjFunction;

class LBFGS {

  public:
    LBFGS(double inTol, int inIter, const std::vector<double> &x, ObjFunction *inObjFunc):
      tol(inTol), iter(inIter), x0(x), objFunc(inObjFunc) 
  {
    //setInitialValue(x0);
  }

    ~LBFGS() {}

  public:
    //void setInitialValue(std::vector<double> x);
    std::vector<double> getCurrentX();
    double getFvalueAfter();
    double lineSearch(std::vector<double> &xold, std::vector<double> &g, std::vector<double> &direction, double stpmax, Mesh *mesh);
    void moveArrayToLeft(std::vector<double> a[], int r);
    bool run(Mesh *mesh);

  public:
    double tol;
    int iter;
    std::vector<double> x0;
    ObjFunction *objFunc;
    std::vector<double> currentX;
    double fValAfter;

  private:
    int r = 50;
}; 

class ObjFunction
{
  public:
    ObjFunction(){};
    //int getSpaceDim() = 0;
    virtual double getTol() = 0;
    double getValue(Mesh *mesh, const std::vector<double> &x) {
      //getValue should get an input of element
      //note that in checkvalidity i am computing it on a all mesh level
      //using a parallel for
      //getvalue should return -ve of worst quality at the cavity level

      return askWorstQuality_2d(mesh, LOs(2,0,1), 2); //2 tri annular mesh test

      /*
      if (dim == 3) calc_crvQuality_3d(mesh);
      // does the following //
      auto rv2v = mesh->ask_down(3, 0).ab2b;
      auto re2e = mesh->ask_down(3, 1).ab2b;
      auto rf2f = mesh->get_adj(3, 2).ab2b;
      auto ev2v = mesh->get_adj(1, 0).ab2b;
      auto vertCtrlPts = mesh->get_ctrlPts(0);
      auto edgeCtrlPts = mesh->get_ctrlPts(1);
      auto faceCtrlPts = mesh->get_ctrlPts(2);
      auto nnew_tets = mesh->nregions();
      auto order = mesh->get_max_order();
      OMEGA_H_CHECK(order == 3);
      auto qs = mesh->ask_qualities();

      Write<Real> Q(nnew_tets, -1e-10);

      auto calc_quality = OMEGA_H_LAMBDA (LO n) {
        Few<Real, 60> tet_pts = collect_tet_pts(order, n, ev2v, rv2v, vertCtrlPts
            , edgeCtrlPts, faceCtrlPts, re2e, rf2f);
        Few<Real, 84> nodes_det = getTetJacDetNodes<84>(3, tet_pts);

        auto const minJ = calcMinJacDet(nodes_det);
        auto const maxJ = calcMaxJacDet(nodes_det);
        if (minJ > 0.) {
          Q[n] = std::pow(std::pow((minJ/maxJ), 1./3.)*qs[n], 1.0);
          if (Q[n] < 0.01) printf("low quality %f for element %d\n", Q[n], n);
        }
        else {
          Q[n] = 0.;
          //printf("tet %d qs %f, minJ %f, maxJ %f, Q %f\n",n,qs[n], minJ, maxJ, Q[n]);
        }
      };
      parallel_for(nnew_tets, std::move(calc_quality));

      return Reals(Q);
      // end of crv quality 3d //////
      */

    }
    /* std::vector<double> getInitialGuess() = 0; */
    /* void setNodes(const vector<double> &x) = 0; */
    /* void restoreInitialNodes() = 0; */ // backtracking if qual. worsens
    // TODO :: can we do this once for all the objective functions?
    std::vector<double> getGrad(Mesh *mesh, const std::vector<double> &_x) {
      double h;
      std::vector<double> x = _x;
      double eps = this->getTol();
      std::vector<double> g;
      for (size_t i = 0; i < x.size(); i++) {
	h = abs(x[i]) > eps ? eps * abs(x[i]) : eps;

	// forward diff
	x[i] += h;
	double ff = this->getValue(mesh, x);
	x[i] -= h;

	// backward diff
	x[i] -= h;
	double fb = this->getValue(mesh, x);
	x[i] += h;

	g.push_back( (ff - fb) / 2./ h );
      }
      return(g);
    }
    ~ObjFunction(){};
};

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

double quadratic_test(const VectorXd& x, VectorXd& grad)
{
  const int n = x.size();
  VectorXd d(n);
  for(int i = 0; i < n; i++)
    d[i] = i;

  double f = (x - d).squaredNorm();
  grad.noalias() = 2.0 * (x - d);
  return f;
}

}

#endif
