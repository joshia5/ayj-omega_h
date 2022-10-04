#ifndef OMEGA_H_IRRULE_HPP
#define OMEGA_H_IRRULE_HPP

#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_few.hpp"

#include <iostream>
#include <fstream>

namespace Omega_h {

OMEGA_H_INLINE Few<Real, 44> TetGaussLobatto(const int order = 4) {
  OMEGA_H_CHECK(order == 4);
  Few<Real, 11*4> ir; //11 points, store x,y,z,weights

  //AddTetPoints4(0, 1./14., 343./45000.);
  { 
    Real a = 1./14.;
    Real b = 1.-3.*a;
    Real weight = 343./45000.;
    ir[0*4+0]=a; ir[0*4+1]=a; ir[0*4+2]=a; ir[0*4+3]=weight;
    ir[1*4+0]=a; ir[1*4+1]=a; ir[1*4+2]=b; ir[1*4+3]=weight;
    ir[2*4+0]=a; ir[2*4+2]=b; ir[2*4+1]=a; ir[2*4+3]=weight;
    ir[3*4+0]=b; ir[3*4+1]=a; ir[3*4+2]=a; ir[3*4+3]=weight;
  }

  //AddTetMidPoint;
  { 
    Real a = 0.25;
    Real weight = -74./5625.;
    ir[4*4+0]=a; ir[4*4+1]=a; ir[4*4+2]=a; ir[4*4+3]=weight;
  }

  //AddTetPoints6(5, 0.10059642383320079500, 28./1125.);
  {
    Real a = 0.10059642383320079500;
    Real b = 0.5-a;
    Real weight = 28./1125.;
    ir[5*4+0]=a; ir[5*4+1]=a; ir[5*4+2]=b; ir[5*4+3]=weight;
    ir[6*4+0]=a; ir[6*4+2]=b; ir[6*4+1]=a; ir[6*4+3]=weight;
    ir[7*4+0]=b; ir[7*4+1]=a; ir[7*4+2]=a; ir[7*4+3]=weight;
    swap2(a,b);
    ir[8*4+0]=a; ir[8*4+1]=a; ir[8*4+2]=b; ir[8*4+3]=weight;
    ir[9*4+0]=a; ir[9*4+2]=b; ir[9*4+1]=a; ir[9*4+3]=weight;
    ir[10*4+0]=b; ir[10*4+1]=a; ir[10*4+2]=a; ir[10*4+3]=weight;
  }
 
  return ir;
}

} //end namespace

#endif
