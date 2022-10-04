#ifndef OMEGA_H_IRRULE_HPP
#define OMEGA_H_IRRULE_HPP

#include "Omega_h_array.hpp"

#include <iostream>
#include <fstream>

namespace Omega_h {
/// Class for integration point with weight
   Real x, y, z, weight;
   LO index;

   void Init(LO const i)
   {
      x = y = z = weight = 0.0;
      index = i;
   }
}

OMEGA_H_DEVICE void TetGaussLobatto(const int order = 4) {

  OMEGA_H_CHECK(order == 4);
  IntegrationRule *ir;
  // 11 points - degree 4 (negative weight)
  ir = new IntegrationRule(11);
  Few<Real, 11*4> //11 points,x,y,z,weights

  void Set(const double x1, const double x2, const double x3, const double w)
   { x = x1; y = x2; z = x3; weight = w; }
  // add the permutations of (a,a,b)
   void AddTetPoints3(const int off, const double a, const double b,
                      const double weight)
   {
      IntPoint(off + 0).Set(a, a, b, weight);
      IntPoint(off + 1).Set(a, b, a, weight);
      IntPoint(off + 2).Set(b, a, a, weight);
   }
  // given a, add the permutations of (a,a,a,b), where 3*a + b = 1
   void AddTetPoints4(const int off, const double a, const double weight)
   {
      IntPoint(off).Set(a, a, a, weight);
      AddTetPoints3(off + 1, a, 1. - 3.*a, weight); }
  // add the permutations of (a,a,b,b), 2*(a + b) = 1
  void AddTetPoints6(const int off, const double a, const double weight)
  {
      const double b = 0.5 - a;
      AddTetPoints3(off,     a, b, weight);
      AddTetPoints3(off + 3, b, a, weight);
  }
  void AddTetMidPoint(const int off, const double weight)
   { IntPoint(off).Set(0.25, 0.25, 0.25, weight); }

  ir->AddTetPoints4(0, 1./14., 343./45000.);
  ir->AddTetMidPoint(4, -74./5625.);
  ir->AddTetPoints6(5, 0.10059642383320079500, 28./1125.);
 
  return ir;
} //end namespace
