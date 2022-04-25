// Import necessary stuff
#include "SOMMG.h"

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Copyright (C) 2013-2019,
-------------------------------------------------------------------------- */

SOMMG::SOMMG(Vec x, PetscInt order) {
    
    //da_mesh=NULL;
    reductionOrder=order;
    
    // Reduction matrix as multivector
    VecDuplicateVecs(x, reductionOrder, &QQ);
    
    // Help vectors and matrices
    VecDuplicateVecs(x, reductionOrder, &Qb);
    
    // Reduced systems
    VecCreateSeq(PETSC_COMM_SELF,reductionOrder,&Pr);
    MatCreateSeqDense(PETSC_COMM_SELF,reductionOrder,reductionOrder,NULL,&Kr);
    MatCreateSeqDense(PETSC_COMM_SELF,reductionOrder,reductionOrder,NULL,&Mr);
    MatCreateSeqDense(PETSC_COMM_SELF,reductionOrder,reductionOrder,NULL,&Cr);
    
}

// Destructor
SOMMG::~SOMMG() {
    
    VecDestroyVecs(reductionOrder, &Qb);
    VecDestroyVecs(reductionOrder, &QQ);
    MatDestroy(&Kr);
    MatDestroy(&Mr);
    MatDestroy(&Cr);
    VecDestroy(&Pr);
    //MatDestroy(&U);

}

PetscErrorCode SOMMG::ConstructBasis(Mat K, Mat M, Mat C, Vec P, Vec N, PetscScalar om, KSP ksp) {
    PetscErrorCode ierr = 0;

    PetscPrintf(PETSC_COMM_WORLD, "Constructing SOMMG reduction basis\n");
    
    // declare and allocate variables (BC vec, moment matching vec, RHS for basis construction)
    Vec NI, momentVec, massVec, fvec;
    VecDuplicate(N, &fvec);
    VecDuplicate(N, &NI);
    VecDuplicate(N, &momentVec);
    VecDuplicate(N, &massVec);
    Mat P0;
    MatDuplicate(K,MAT_COPY_VALUES,&P0);
    MatAXPY(P0,-om*om,M,SAME_NONZERO_PATTERN); // quasi static: 0
    
    // Container for the moment matching
    PetscScalar *Utemp = new PetscScalar[reductionOrder*reductionOrder];
    
    // Zero all containers
    for (PetscInt i=0;i<reductionOrder;i++){
        VecSet(QQ[i], 0.0); 
        VecSet(Qb[i], 0.0);
        for (PetscInt j=0;j<reductionOrder;j++){
            Utemp[i*reductionOrder+j] = 0.0;
        }
    }
    
    // Help vector for BCs on full problem
    VecSet(NI, 1.0);
    VecAXPY(NI, -1.0, N);
    
    // start time
    double t1 = MPI_Wtime();
    
    // Set BCs on full matrices - and add ones to diag of K
    MatDiagonalScale(P0, N, N);
    MatDiagonalScale(M, N, N);
    MatDiagonalScale(C, N, N);
    MatDiagonalSet(P0, NI, ADD_VALUES);

    // Zero out possible loads in the P that coincide with Dirichlet conditions
    VecPointwiseMult(P, P, N);
    
    // STEP 0: Create the first reduction vector, QQ[0[]
    // Solve at the expansion frequency
    ierr = KSPSetOperators(ksp, P0, P0); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);
    ierr = KSPSolve(ksp, P, Qb[0]);
    
    // Get solver information and print to screen
    PetscInt    niter;
    PetscScalar rnorm;
    KSPGetIterationNumber(ksp, &niter);
    KSPGetResidualNorm(ksp, &rnorm);
    PetscReal RHSnorm;
    ierr = VecNorm(P, NORM_2, &RHSnorm); CHKERRQ(ierr);
    rnorm = rnorm / RHSnorm;
    double t2 = MPI_Wtime();
    PetscPrintf(PETSC_COMM_WORLD, "SOMMG Reduction basis %i of %i  solve:  iter: %i, rerr.: %e, time: %f\n",1,reductionOrder, niter, rnorm, t2 - t1);
    
    // Start building U-matrix, i.e. U(1,1) = norm(Qb(:,1));
    ierr = VecNorm(Qb[0], NORM_2, &Utemp[0]); CHKERRQ(ierr);
    
    // Q(:,1) = Qb(:,1)/U(1,1);
    VecAXPY(QQ[0],1.0/Utemp[0],Qb[0]);
        
    // LOOP TO CONSTRUCT THE REST OF THE BASIS
    for (PetscInt j=1;j<reductionOrder;j++){
        
       // PetscPrintf(PETSC_COMM_WORLD, "Start iter: %i of %i \n",j, reductionOrder-1);
        
        // STEP 1: Create the RHS for the next basis vector
        // fvec = -1.0*P1*Q[j-1]; Note that this is the same for all  P1 = -2*om*M
        MatMult(M,QQ[j-1],fvec);
        VecScale(fvec,(-1.0)*(-2.0*om));
                
        // STEP 2: If j>1, introduce the moment matching correction to the RHS
        if (j>1) {
            // Compute RHS with moment matching correction, i.e.
            // fvec = -P0 * Q[j-1] - P2 * ( Q[:,0:j-2] * ( U[1:j-1,1:j-1]^-1 * eUnit ) );
            // Note: First part already computed - find the rest and subtract from fvec
            // Second part in steps - from right to left:
        
            // Declare and allocate data containers (all local and small!)
            Vec eUnit, bvec;
            Mat U; 
            VecCreateSeq(PETSC_COMM_SELF,j-1,&eUnit);
            VecCreateSeq(PETSC_COMM_SELF,j-1,&bvec);
            MatCreateSeqDense(PETSC_COMM_SELF,j-1,j-1,NULL,&U);
            
            // MOMENT STEP 1: Create the e[j-2] unit vector
            VecSet(eUnit,0.0);
            VecSetValue(eUnit,j-2,1.0,INSERT_VALUES);
            VecAssemblyBegin(eUnit);
            VecAssemblyEnd(eUnit);
            
            // MOMENT STEP 2: insert Utemp[1:j-1,1:j-1] into U
            MatZeroEntries(U);
            for (PetscInt k=1;k<j;k++){
                for (PetscInt h=1;h<j;h++){
                    MatSetValue(U,k-1,h-1,Utemp[k*reductionOrder+h],INSERT_VALUES);
                }
            }
            MatAssemblyBegin(U,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(U,MAT_FINAL_ASSEMBLY); 
           
            // MOMENT STEP 3: Solve for bvec = U[1:j-1,1:j-1]^-1 * eUnit
            t1 = MPI_Wtime();
            KSP ksploc;
            KSPCreate(PETSC_COMM_SELF, &(ksploc));
            KSPSetFromOptions(ksploc); // Note - all runtime settings overwritten next!
            KSPSetType(ksploc, KSPPREONLY);
            PC pcloc;
            KSPGetPC(ksploc, &pcloc);
            PCSetFromOptions(pcloc); // Note - all runtime settings overwritten next!
            PCSetType(pcloc, PCLU);
            KSPSetOperators(ksploc, U, U);
            KSPSetUp(ksploc);                    
            KSPSolve(ksploc, eUnit, bvec);
            t2 = MPI_Wtime();
         
            // MOMENT STEP 4: Multiply and get momentVec = Q[:,0:j-2] * bvec
            // extract values from bvec
            PetscScalar *pbvec;
            VecGetArray(bvec,&pbvec);
            // Zero the momentVec
            VecSet(momentVec, 0.0);
            for (PetscInt k=0;k<j-1;k++){
                VecAXPY(momentVec,pbvec[k],QQ[k]);
            }
            VecRestoreArray(bvec,&pbvec);
            
            // MOMENT STEP 5: multiply on mass matrix
            MatMult(M,momentVec,massVec);
            
            // MOMENT STEP 6: Subtract from fvec
            VecAXPY(fvec,(-1.0)*(-1,0),massVec);
            
            // Clean up
            KSPDestroy(&ksploc);
            MatDestroy(&U);
            VecDestroy(&eUnit);
            VecDestroy(&bvec);
        
        }

        // STEP 3: Solve
        // Apply boundary conditions, i.e. fvec = Null*fvec;
        VecPointwiseMult(fvec, fvec, N);

        // Solve for fvec, i.e. Qb(j) = dKeff\fvec;
        t1 = MPI_Wtime();
        ierr = KSPSolve(ksp, fvec, Qb[j]);
        t2 = MPI_Wtime();
        KSPGetIterationNumber(ksp, &niter);
        KSPGetResidualNorm(ksp, &rnorm);
        ierr = VecNorm(fvec, NORM_2, &RHSnorm); CHKERRQ(ierr);
        rnorm = rnorm / RHSnorm;
        PetscPrintf(PETSC_COMM_WORLD, "SOMMG Reduction basis %i of %i  solve:  iter: %i, rerr.: %e, time: %f\n", j+1,reductionOrder,niter, rnorm, t2 - t1);
        
        // STEP 4: Orthonormalize 
        for (PetscInt i=0;i<j;i++){
            // Utemp[i,j] = QQ_i^T*Qb_j;
            VecDot(QQ[i],Qb[j],&( Utemp[i*reductionOrder+j] ) );
            //PetscPrintf(PETSC_COMM_WORLD, "Index %i, Unorm %e \n",j,Utemp[i*reductionOrder+j]);
            // Qb_j = Qb_j - Utemp[i,j]*QQ_i;
            VecAXPY(Qb[j],-1.0*Utemp[i*reductionOrder+j],QQ[i]);
        }
       
        // STEP 5: Update diagonal of moment matching matrix and compute the reduction basis
        // Utemp[j,j] = norm2(Qb_j);
        ierr = VecNorm(Qb[j], NORM_2, &(Utemp[j*reductionOrder+j])); CHKERRQ(ierr);
        // QQ_j = Qb_j / Utemp[j,j];
        VecAXPY(QQ[j],1.0/Utemp[j*reductionOrder+j],Qb[j]);
        
    } // End basis construction loop
    
    // Clean up
    VecDestroy(&fvec);
    VecDestroy(&NI);
    VecDestroy(&momentVec);
    VecDestroy(&massVec);
    MatDestroy(&P0);
    delete [] Utemp;

    // Put below part in its own method - can be used by all Krylov reduction methods
    {
        // Create reduced matrices and RHS
        MatZeroEntries(Kr);
        MatZeroEntries(Mr);
        MatZeroEntries(Cr);
        PetscScalar dotval;
        Vec tmp;
        VecDuplicate(P,&tmp);
        for (PetscInt i=0;i<reductionOrder;i++){
            VecDot(P,QQ[i],&dotval);
            VecSetValues(Pr,1,&i,&dotval,INSERT_VALUES);
            for (PetscInt j=0;j<reductionOrder;j++){
                // Entry for K
                MatMult(K,QQ[j],tmp);
                VecDot(tmp,QQ[i],&dotval);
                MatSetValues(Kr,1,&i,1,&j,&dotval,INSERT_VALUES);
                // Entry for M
                MatMult(M,QQ[j],tmp);
                VecDot(tmp,QQ[i],&dotval);
                MatSetValues(Mr,1,&i,1,&j,&dotval,INSERT_VALUES);
                // Entry for C
                MatMult(C,QQ[j],tmp);
                VecDot(tmp,QQ[i],&dotval);
                MatSetValues(Cr,1,&i,1,&j,&dotval,INSERT_VALUES);
            }
        }
        MatAssemblyBegin(Kr,MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Mr,MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(Cr,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Kr,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Mr,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Cr,MAT_FINAL_ASSEMBLY);
        
        VecDestroy(&tmp);
    }
 

   
    return ierr;
}
