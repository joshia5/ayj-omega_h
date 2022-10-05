#ifndef OMEGA_H_IRRULE_HPP
#define OMEGA_H_IRRULE_HPP

#include "Omega_h_array.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_few.hpp"

#include <iostream>
#include <fstream>

namespace Omega_h {


  const DenseMatrix &Wideal =
    Geometries.GetGeomToPerfGeomJac(fe.GetGeomType());
  MFEM_ASSERT(Wideal.Height() == Jtr.SizeI(), "");
  MFEM_ASSERT(Wideal.Width() == Jtr.SizeJ(), "");

        for (int i = 0; i < ir.GetNPoints(); i++) { Jtr(i) = Wideal; }
        break;
  Few<Real, 9999> Jtr; //

    return Jtr;

} //end namespace

#endif
