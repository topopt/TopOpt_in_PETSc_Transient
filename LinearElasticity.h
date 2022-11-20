#ifndef __LINEARELASTICITY__
#define __LINEARELASTICITY__

#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc.h>
#include <petsc/private/dmdaimpl.h>
#include <MPIIO.h>
#include <SOMMG.h>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
 Updated: June 2019, Niels Aage
 Copyright (C) 2013-2019,

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

class LinearElasticity {

  public:

    // Constructor
    //LinearElasticity(DM da_nodes, DM da_elem);
    LinearElasticity(DM da_nodes, DM da_elem,PetscScalar nu_in, PetscInt nlvls_in, PetscInt ncp_in, PetscInt tsc_in, PetscScalar t1_in, PetscInt ti_in,PetscInt objective_in, PetscInt stiffnessInterpolation_in,PetscInt massInterpolation_in, PetscScalar t0_in,PetscInt boundaries_in,PetscInt loadloc_in,PetscBool reduction_in,PetscInt basis_in, PetscScalar dalpha_in, PetscScalar dbeta_in,PetscScalar omega0_in,PetscInt optimizationRegion, PetscScalar V0); 
	 
	 // Destructor
    ~LinearElasticity();

    //  Compute objective and constraints and sensitivities at once: GOOD FOR
    //  SELF_ADJOINT PROBLEMS
    //PetscErrorCode ComputeObjectiveConstraintsSensitivities(PetscScalar* fx, PetscScalar* gx, Vec dfdx, Vec dgdx, Vec xPhys, PetscScalar Emin,        PetscScalar Emax,PetscScalar Rmin, PetscScalar Rmax,PetscScalar penalK,PetscScalar penalM, PetscScalar volfrac, PetscScalar dalpha, PetscScalar dbeta);

//    PetscErrorCode ComputeObjective(PetscScalar* fx, Vec xPhys ,PetscScalar Emin,   PetscScalar Emax,PetscScalar Rmin, PetscScalar Rmax,PetscScalar penalK,PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta);
//    PetscErrorCode ComputeSensitivities(Vec dfdx, Vec xPhys, PetscScalar Emin,        PetscScalar Emax,PetscScalar Rmin, PetscScalar Rmax,PetscScalar penalK,PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta);
    //PetscErrorCode ComputeObjectiveSensitivities(PetscScalar* fx, Vec dfdx, Vec xPhys, PetscScalar Emin,        PetscScalar Emax,PetscScalar Rmin, PetscScalar Rmax,PetscScalar penalK,PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta);

PetscErrorCode ComputeObjectiveSensitivities(PetscScalar* fx,Vec dfdx,Vec xPhys,  PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta,PetscBool write,MPIIO * output);

    // Solve the FE problem
    PetscErrorCode SolveState(Vec xPhys, PetscScalar Emin, PetscScalar Emax,PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta);


    // Solve the lagrange multipliers
//    PetscErrorCode SolveLagrangen(PetscInt c);
		PetscErrorCode SolveLagrangen(PetscInt c,Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta);
//    PetscErrorCode SolveLagrange0(PetscInt c);
PetscErrorCode SolveLagrange0(PetscInt c,Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta);

    // Compute objective and constraints for the optimiation
    //PetscErrorCode ComputeObjectiveConstraints(PetscScalar* fx, PetscScalar* gx, Vec xPhys, PetscScalar Emin,PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK,PetscScalar penalM, PetscScalar volfrac, PetscScalar dalpha, PetscScalar dbeta);

    // Compute sensitivities
    //PetscErrorCode ComputeSensitivities(Vec dfdx, Vec dgdx, Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar penal,PetscScalar volfrac); // needs ....


    // Restart writer
    PetscErrorCode WriteRestartFiles();

    // Get pointer to the FE solution
    // Vec GetStateFieldU() { return (&U); };
    // Vec GetStateFieldDU() { return (&DU); };
    // Vec GetStateFieldDDU() { return (&DDU); };
    Vec GetStateFieldU(PetscInt j) { return (U[j]); };
    Vec GetStateFieldDU(PetscInt j) { return (DU[j]); };
    Vec GetStateFieldDDU(PetscInt j) { return (DDU[j]); };

    // Get pointer to DMDA
    DM GetDM() { return (da_nodal); };

    // Logical mesh
    DM da_nodal; // Nodal mesh

    PetscInt Getts(){ return (Ts); }

    // Sets up solvers, completes initial conditions DDU0 and sets system matrix for continued time integration
    PetscErrorCode CompleteInitialState(Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta,PetscBool write,MPIIO *output2);

    PetscErrorCode GenerateCheckPoints(PetscBool write, Vec xPhys,MPIIO *output2);

    PetscErrorCode VolumeConstraint(PetscScalar* gx,Vec dgdx, Vec xPhys,PetscScalar volfrac);


  private:
    // Working dir
    std::string wd0;
	
		 //SOAR* reductionmodel;
		 SOMMG* reductionmodel;

    // Logical mesh
    PetscInt    nn[3]; // Number of nodes in each direction
    PetscInt    ne[3]; // Number of elements in each direction
    PetscScalar xc[6]; // Domain coordinates
	 const PetscScalar pi = 3.14;

    // Linear algebra
	 PetscInt	 objective; // choose between multiple objective functions 
	 PetscScalar dalpha, dbeta; 
	 PetscInt basis ;
	 PetscBool  reduction;
	 PetscScalar omega0;
	 PetscInt    boundaries; //Boundary conditions selector 
	 PetscInt    optimizationRegion;
    PetscInt    loadloc; // Load position
	 Mat         M_full,C_full,K_full;           // Global stiffness matrix
	 Mat         M,C,K;           // Global stiffness matrix
    Mat         A,B1,B2,B3,dRdS;        // Jacobian matrix
    Vec         U0,U0_full,DU0,DU0_full,DDU0,DDU0_full;   // Displacement vector
    Vec         L0,L0_full,DL0,DL0_full,DDL0,DDL0_full, LT;           // Displacement vector
    Vec         UM1, DUM1, DDUM1; // "n minus one " state vectores
    Vec         R;          // residual vector
    Vec         Noriginal,N;          // Dirichlet vector (used when imposing BCs)
    Vec         P;          // Elemet identificator vector 
    Vec         *U,*DU,*DDU;           // Displacement vector
    Vec         *U_full,*DU_full,*DDU_full;           // Displacement vector
    Vec         *L,*DL,*DDL;           // Lagrange arrrays
    Vec         *L_full,*DL_full,*DDL_full;           // Lagrange arrrays
	 Vec			 HDU,HDU_full,HDDU,HDDU_full,HDDDU,HDDDU_full;	   // lagrange solve loads 
    Vec         RHSLocation,RHS0,RHS;       // Load vector
    Vec         Lp1,DLp1,DDLp1; // Lagrange vector saves
    PetscScalar hdu0, hddu0, hdddu0; 
	 PetscScalar ME[24 * 24]; // Element stiffness matrix
    PetscScalar CE[24 * 24]; // Element stiffness matrix
    PetscScalar KE[24 * 24]; // Element stiffness matrix
    PetscScalar a[6];
	 PetscScalar s, ds; // Elementwise element scalings 
	 PetscScalar h0, dfx0; // elementwise contribution to the objective function

    // Timeintegration
    PetscInt Ts; // NUmber of timesteps
    PetscInt ncp;
    PetscInt tsc;
    PetscScalar dt;
    PetscReal t0; //
    PetscReal t1; //
    PetscInt ti; //
    PetscInt t; //
    PetscScalar epsi;

	// interpolation properties
	PetscScalar mfac;
	PetscScalar dmfac;
	PetscScalar kfac;
	PetscScalar dkfac;
	PetscInt stiffnessInterpolation;
	PetscInt massInterpolation;

    // Solver
    KSP         ksp,ksp1,ksp2; // Pointer to the KSP object i.e. the linear solver+prec
    PetscInt    nlvls;
    PetscScalar nu; // Possions ratio
    PetscScalar V0; // Possions ratio

    // Set up the FE mesh and data structures
    PetscErrorCode SetUpBC(DM da_nodes,PetscInt boundaries);
    PetscErrorCode SetUpTimeIntegration(PetscInt ti);
    PetscErrorCode SetUpOptimizationRegion(DM da_nodes, DM da_elem);

    // Assemble the stiffness matrix
    PetscErrorCode AssembleStiffnessMatrix(Vec xPhys, PetscScalar min, PetscScalar max, PetscScalar penal);
    PetscErrorCode AssembleMassMatrix(Vec xPhys, PetscScalar min, PetscScalar max, PetscScalar penal);
    PetscErrorCode AssembleDampingMatrix(PetscScalar dalpha, PetscScalar dbeta);

	// Interpolation
	PetscErrorCode InterpolateStiffness(PetscScalar Emax, PetscScalar Emin,PetscScalar penalK,PetscScalar xp);
	PetscErrorCode InterpolateMass(PetscScalar Rmax, PetscScalar Rmin,PetscScalar penalM, PetscScalar xp);


// Put the objective function evaluation in this function 
PetscErrorCode EvaluateObjective(PetscScalar kfac, PetscScalar uKu, PetscScalar mfac,PetscScalar duMdu,PetscScalar u,PetscScalar du, PetscScalar uu, PetscScalar dudu,PetscScalar ddu, PetscScalar dduddu,PetscScalar duCdu);
PetscErrorCode EvaluateObjectiveGradient( const  PetscScalar dkfac, const PetscScalar uKu, const PetscScalar dmfac, const PetscScalar duMdu, const PetscScalar duCdu, const PetscScalar DxduCdu);
//PetscErrorCode GradObjeciveWrtState(const PetscScalar dt, const PetscScalar u,  const PetscScalar du, const PetscScalar ddu, const PetscScalar K, const PetscScalar M);
PetscErrorCode GradObjectiveWrtState(PetscInt k, PetscScalar dt, PetscScalar *up,PetscScalar *dup,PetscScalar *ddup,PetscInt *edof,PetscScalar kfac,PetscScalar mfac,PetscScalar dalpha,PetscScalar dbeta);

PetscErrorCode EvaluateObjectiveScale(PetscScalar xp);

// Time integrate checkpoint block
    PetscErrorCode SolveStateBlock(Vec U0, Vec DU0, Vec DDU0,PetscInt c);

    // Load functions
    PetscErrorCode Load(PetscScalar t);
    PetscErrorCode LoadLocation();
	 PetscErrorCode LoadPoint(PetscInt i, PetscScalar xcordp,PetscScalar ycordp,PetscScalar zcordp,PetscInt axis, PetscScalar t);

    // Start the solver
    PetscErrorCode SetUpSolver1();
    PetscErrorCode SetUpSolver2();
    PetscErrorCode LoadSolutionAsU0();

    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]);

    // Methods used to assemble the element stiffness matrix
    PetscInt    Hex8IsoparametricKE(PetscScalar* X, PetscScalar* Y, PetscScalar* Z, PetscScalar nu, PetscInt redInt,PetscScalar* ke);
    PetscInt    Hex8IsoparametricME(PetscScalar* X, PetscScalar* Y, PetscScalar* Z, PetscInt redInt,PetscScalar* me);
    PetscInt Hex8IsoparametricCE(PetscScalar dalpha, PetscScalar dbeta,PetscScalar* ce);

    PetscScalar Dot(PetscScalar* v1, PetscScalar* v2, PetscInt l);
    void        ShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta, PetscScalar* S);
    void        ShapeFunctionsVector(PetscScalar* S,PetscScalar* vS);
    void        DifferentiatedShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta, PetscScalar* dNdxi,PetscScalar* dNdeta, PetscScalar* dNdzeta);
    PetscScalar Inverse3M(PetscScalar J[][3], PetscScalar invJ[][3]);

//PetscErrorCode ComputeObjectiveSensitivitiesElementLoopn(PetscScalar* fx, Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc, Vec Lloc, Vec Ubloc, Vec DUbloc,Vec DDUbloc,Vec N);

//PetscErrorCode ComputeObjectiveSensitivitiesElementLoop0(PetscScalar* fx, Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc, Vec DDLloc, Vec N);



PetscErrorCode ComputeSensitivitiesElementLoopn(Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc, Vec Lloc, Vec Ubloc, Vec DUbloc,Vec DDUbloc,Vec N,Vec P);

PetscErrorCode ComputeSensitivitiesElementLoop0(Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc, Vec DDLloc, Vec N,Vec P);


PetscErrorCode ComputitiesElementLoopn(Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc, Vec Lloc, Vec Ubloc, Vec DUbloc,Vec DDUbloc,Vec N);


PetscErrorCode ComputeObjectiveElementLoop(PetscScalar* fx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc, Vec DDUloc );

PetscErrorCode ComputeGradientObjectiveElementLoop( Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec U, Vec DU,Vec DDU);


    // Restart
    PetscBool   restart, flip;
    std::string filename00, filename01;

    // File existence
    inline PetscBool fexists(const std::string& filename) {
        std::ifstream ifile(filename.c_str());
        if (ifile) {
            return PETSC_TRUE;
        }
        return PETSC_FALSE;
    }
};

#endif
