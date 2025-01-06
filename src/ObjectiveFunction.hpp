#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <iostream>
#include <vector>

//#include "Omega_h_mesh.hpp"
#include "Omega_h_beziers.hpp"

#include <Eigen/Core>
#include <LBFGS.h>
#include <iostream>
//namespace Omega_h {}

class ObjectiveFunction
{
  private:
    Mesh* mesh;
    int n;
  public:
    ObjectiveFunction(int n_, Mesh* mesh_) : mesh(mesh_) {}
    double operator() (const Eigen::VectorXd& x, Eigen::VectorXd& grad)
    {
      //forget about Mesh for now, just do it for a given set of input
      //doubles, look at the body of the askWorstQual and come in at the point
      //where an array or whatever containing coordinate locations of control
      //points is provided to the detJ calc function

      grad[1] = 20;
      grad[0]     = -2.0;

      return askWorstQuality_2d(mesh, LOs(mesh->nelems(),0,1), 2);

    }

    /*
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
    */
    ~ObjectiveFunction(){};
};

#endif
