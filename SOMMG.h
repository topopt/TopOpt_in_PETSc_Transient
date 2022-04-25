#ifndef SOMMG_H
#define SOMMG_H

// Include necessary libraries
#include <petsc.h>
#include <iostream>
#include <mpi.h>

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Copyright (C) 2020
-------------------------------------------------------------------------- */

class SOMMG {
  public:
      
    // ------------- CONSTRUCT/DESTRUCT -------------------------------
     SOMMG(Vec x, PetscInt order);
    ~SOMMG();

    // ------------- REDUCED ORDER MATRIX -----------------------------
    Vec *QQ, Pr;                              // Reduction "matrix" QQ and reduced load vector Pr
    Mat Mr,Kr,Cr;                             // Reduced matrices
    
    // ------------- METHODS ------------------------------------------
    PetscErrorCode ConstructBasis(Mat K, Mat M, Mat C, Vec P, Vec N, PetscScalar om, KSP ksp);
 
  private:
    
    PetscInt reductionOrder;                  // degree of expansion
    Vec *Qb;                                 // Orthonormalization "matrix"

};

#endif
