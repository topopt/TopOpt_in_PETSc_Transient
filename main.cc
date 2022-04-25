#include <Filter.h>
#include <LinearElasticity.h>
#include <MMA.h>
#include <MPIIO.h>
#include <TopOpt.h>
#include <mpi.h>
#include <petsc.h>
/*
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

Updated: June 2019, Niels Aage
Copyright (C) 2013-2019,

Disclaimer:
The authors reserves all rights but does not guaranty that the code is
free from errors. Furthermore, we shall not be liable in any event
caused by the use of the program.
*/

static char help[] = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";


int main(int argc, char* argv[]) {

    // Error code for debugging
    PetscErrorCode ierr;

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
    TopOpt* opt = new TopOpt();

    // // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes, opt->da_elem,opt->nu,opt->nlvls,opt->ncp,opt->tsc,opt->t1,opt->ti,opt->objective,opt->sInt,opt->mInt,opt->t0,opt->boundaries,opt->loadloc,opt->reduction,opt->basis,opt->dalpha,opt->dbeta,opt->omega0,opt->optimizationRegion,opt->V0); 

    // STEP 4: VISUALIZATION USING VTK
    MPIIO* output1 = NULL;
	MPIIO* output2 = NULL;
	MPIIO* output3 = NULL;
    if (opt->outputDesign) {
	  output1 = new MPIIO(opt->da_nodes, 0, "", 3, "x, xTilde, xPhys","outputDesign_00000.dat");
	 }
	 if (opt->outputPhysics0) {
    output2 = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 1, "xPhys","outputPhysics0_00000.dat");
	 }
	 if (opt->outputPhysics1){
    output3 = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 1, "xPhys","outputPhysics1_00000.dat");
	 }

	 // ALLOCATE FILTER
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);
	
    // STEP 5: THE OPTIMIZER MMA
    MMA*     mma;
    PetscInt itr = 0;
    opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !
	 mma->SetAsymptotes(0.1, 0.6, 1.05);

    // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
    filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);

	 // STEP 7: OPTIMIZATION LOOP
    PetscScalar ch = 1.0;
    double      time1, t2 ;
    PetscBool  firstIterWrite;
    PetscScalar chlim = 0.01;

	//Indicate if this is the first iteration 
	firstIterWrite = PETSC_TRUE;
   

	 while (itr < opt->maxItr && ch > chlim) {
		
	     // Update iteration counter
        itr++;

        // start timer
        time1 = MPI_Wtime();

        // computes the initial accelerations, sets up K, M and C, and sets the KSP-operator for time integration
        physics->CompleteInitialState(opt->xPhys, opt->Emin, opt->Emax,opt->Rmin,opt->Rmax, opt->penalK,opt->penalM,opt->dalpha,opt->dbeta,firstIterWrite,output2);

        // Perform timeintegration and write checkpoints. Writes to disk.
		  if (opt->ncp>1) {
		   	physics->GenerateCheckPoints(firstIterWrite,opt->xPhys,output2);
		  }

		  // Compute objective function and sensitivities 
        physics->ComputeObjectiveSensitivities(&(opt->fx),opt->dfdx, opt->xPhys, opt->Emin, opt->Emax,opt->Rmin,opt->Rmax, opt->penalK,opt->penalM,opt->dalpha,opt->dbeta,firstIterWrite,output2);    //

        // Compute volume constraints
        physics->VolumeConstraint(&(opt->gx[0]),opt->dgdx[0], opt->xPhys,opt->volfrac);
        
        // Filter sensitivities (chainrule)
        filter->Gradients(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta, opt->eta);
        
		//VecView(opt->dfdx,PETSC_VIEWER_STDOUT_WORLD);
		  // Make sure nothing no visual outputs is written in opt.loop from this point
        firstIterWrite = PETSC_FALSE; 

        // design output
        if (itr % 1 == 0 && output1!=NULL) {
            output1->WriteVTK(physics->da_nodal,opt->x, opt->xTilde, opt->xPhys, itr);
        }

        // Discreteness measure
        PetscScalar mnd = filter->GetMND(opt->xPhys);
        
		  // Compute objective scale
        if (itr == 1) {
            opt->fscale = 1.0 / opt->fx;
        }

        // Scale objectie and sens
        opt->fx = opt->fx * opt->fscale;
        VecScale(opt->dfdx, opt->fscale);

        // Sets outer movelimits on design variables
        mma->SetOuterMovelimit(opt->Xmin, opt->Xmax, opt->movlim, opt->x, opt->xmin, opt->xmax);

        // Update design by MMA
		  if (opt->minmax==1) { // try to maximize instead of minimize objective 
			VecScale(opt->dfdx,-1);
			}
        mma->Update(opt->x, opt->dfdx, opt->gx, opt->dgdx, opt->xmin, opt->xmax);

		 // Inf norm on the design change
        ch = mma->DesignChange(opt->x, opt->xold);

        // stop timer
        t2 = MPI_Wtime();

		  // Print to screen
        PetscPrintf(PETSC_COMM_WORLD,"# It.: %03i, True fx: %4.10e, Scaled fx: %4.2f, gx[0]: %+02.2f, ch: %02.2e, mnd: %2.2f, time: %02.2f\n", itr, opt->fx / opt->fscale, opt->fx, opt->gx[0], ch, mnd, t2 - time1);
        
		  // Increase beta if needed
		  PetscBool changeBeta = PETSC_FALSE;
        if (opt->projectionFilter && (itr >= 200 || ch < 0.01)) {
			 
					changeBeta = filter->IncreaseBeta(&(opt->beta), opt->betaFinal, opt->gx[0], itr, ch);
					if (changeBeta  ){
						ch = 1.0; // Ensure the optimization does not get stuck unintentionally. 
						}
			
		  }

		  // ensure that optimization is not terminated if constraints are violated 
		  if (opt->gx[0]>0){
		     ch = 1;
		  }

        // Filter design field
        filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);

       // Dump data needed for restarting code at termination
        if (itr % 10  == 0 || changeBeta) {
            opt->WriteRestartFiles(&itr, mma);
            physics->WriteRestartFiles();
        }
    }

    //Write restart WriteRestartFiles
    opt->WriteRestartFiles(&itr, mma); 
    physics->WriteRestartFiles(); 

    //PetscPrintf(PETSC_COMM_WORLD,"Complete_final_state...\n");
    physics->CompleteInitialState(opt->xPhys, opt->Emin, opt->Emax,opt->Rmin,opt->Rmax, opt->penalK,opt->penalM,opt->dalpha,opt->dbeta,PETSC_TRUE,output3);
    
    // Perform timeintegration and write checkpoints. Writes to disk.
	if (opt->ncp>1) {
			physics->GenerateCheckPoints(PETSC_TRUE,opt->xPhys,output3);
	}

		// Compute objective function and sensitivities 
	physics->ComputeObjectiveSensitivities(&(opt->fx),opt->dfdx, opt->xPhys, opt->Emin, opt->Emax,opt->Rmin,opt->Rmax, opt->penalK,opt->penalM,opt->dalpha,opt->dbeta,PETSC_TRUE,output3);    //

    // // STEP 7: CLEAN UP AFTER YOURSELF
    delete mma;
    if (opt->outputDesign){
        delete output1;
    }
	if (opt->outputPhysics0){
        delete output2;
    }
    if (opt->outputPhysics1){
        delete output3;
    }
    delete filter;
    delete physics;
    delete opt;

	// Termination message
	PetscPrintf(PETSC_COMM_WORLD,"# Optimization has terminated\n");

    // Finalize PETSc / MPI
    PetscFinalize();
    return 0;
}
