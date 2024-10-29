#ifndef LBFGS_H
#define LBFGS_H

#include <iostream>
//#include <pcu_util.h>
#include <vector>

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
    double lineSearch(std::vector<double> &xold, std::vector<double> &g, std::vector<double> &direction, double stpmax);
    void moveArrayToLeft(std::vector<double> a[], int r);
    bool run();

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
    virtual ~ObjFunction(){};
    virtual int getSpaceDim() = 0;
    virtual double getTol() = 0;
    virtual double getValue(const std::vector<double> &x) = 0;
    /* virtual std::vector<double> getInitialGuess() = 0; */
    /* virtual void setNodes(const vector<double> &x) = 0; */
    /* virtual void restoreInitialNodes() = 0; */
    // TODO :: can we do this once for all the objective functions?
    std::vector<double> getGrad(const std::vector<double> &_x)
    {
      double h;
      std::vector<double> x = _x;
      double eps = this->getTol();
      std::vector<double> g;
      for (size_t i = 0; i < x.size(); i++) {
	h = abs(x[i]) > eps ? eps * abs(x[i]) : eps;

	// forward diff
	x[i] += h;
	double ff = this->getValue(x);
	x[i] -= h;

	// backward diff
	x[i] -= h;
	double fb = this->getValue(x);
	x[i] += h;

	g.push_back( (ff - fb) / 2./ h );
      }
      return(g);
    }
};

}

#endif
