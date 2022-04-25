#include <LinearElasticity.h>
#include <math.h>
#include <SOMMG.h>
/*
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

Disclaimer:
The authors reserves all rights but does not guaranty that the code is
free from errors. Furthermore, we shall not be liable in any event
caused by the use of the program.
*/

//LinearElasticity::LinearElasticity(DM da_nodes, DM da_elem,PetscScalar nu_in,PetscInt stiffnessInterpolation_in,PetscInt massInterpolation_in) {
 LinearElasticity::LinearElasticity(DM da_nodes, DM da_elem,PetscScalar nu_in, PetscInt nlvls_in, PetscInt ncp_in, PetscInt tsc_in, PetscScalar t1_in, PetscInt ti_in,PetscInt objective_in, PetscInt stiffnessInterpolation_in,PetscInt massInterpolation_in, PetscScalar t0_in, PetscInt boundaries_in,PetscInt loadloc_in,PetscBool reduction_in, PetscInt basis_in,PetscScalar dalpha_in, PetscScalar dbeta_in,PetscScalar omega0_in,PetscInt optimizationRegion_in, PetscScalar V0_in) { 
//LinearElasticity::LinearElasticity(DM da_nodes, DM da_elem) {
    // Set pointers to null
    M_full   = NULL;
    C_full   = NULL;
    K_full   = NULL;
    M   = NULL;
    C   = NULL;
    K   = NULL;
    U0   = NULL;
    U0_full   = NULL;
    DU0_full   = NULL;
    DDU0_full   = NULL;
    DU0   = NULL;
    DDU0   = NULL;
    L0   = NULL;
    L0_full   = NULL;
	 LT = NULL;
    DL0   = NULL;
    DL0_full   = NULL;
    DDL0   = NULL;
    DDL0_full   = NULL;
    U   = NULL;
    U_full   = NULL;
    L   = NULL;
    L_full   = NULL;
    P   = NULL;
    DU  = NULL;
    DU_full  = NULL;
    DL  = NULL;
    DL_full  = NULL;
    DDU = NULL;
    DDU_full = NULL;
    DDL = NULL;
    DDL_full = NULL;
    Lp1 = NULL;
    DLp1 = NULL;
    DDLp1 = NULL;
    RHS = NULL;
    RHS0 = NULL;
    RHSLocation = NULL;
    HDU = NULL;
    HDDU = NULL;
    HDDDU = NULL;
    HDU_full = NULL;
    HDDU_full = NULL;
    HDDDU_full = NULL;
    N   = NULL;
    Noriginal   = NULL;
    R   = NULL;
    ksp1 = NULL;
    ksp2 = NULL;
    A = NULL;
    dRdS = NULL;
    B1= NULL;
    B2 = NULL;
    B3 = NULL;
	reductionmodel = NULL;

    // Parameters - to be changed on read of variables
	 nu = nu_in;
	 V0 = V0_in;
    nlvls = nlvls_in;//4;  // NUmber of multigrid levels
    ncp = ncp_in;////2; // NUmber of checkpoints
    tsc = tsc_in;//25; // number of timesteps per checkpoint
    t1 = t1_in;//50;   // Final time
    t0 = t0_in;    // Initial time
    ti = ti_in;//1; // must be an integer
	 reduction = reduction_in;
	 basis = basis_in;
	 objective = objective_in;//0;;
	 dalpha = dalpha_in;
	 omega0 = omega0_in;
	 dbeta = dbeta_in;
	 boundaries = boundaries_in;
	 loadloc = loadloc_in;
	 optimizationRegion = optimizationRegion_in;
	 stiffnessInterpolation = stiffnessInterpolation_in;
	 massInterpolation = massInterpolation_in;
    PetscBool flg;

	 // Total number if timesteps 	
	 Ts = ncp*tsc+1;                    
	 
	 // time step size 
    dt = t1/(Ts-1);

    // Identify workdir
    char      wdgrab[PETSC_MAX_PATH_LEN];
    PetscOptionsGetString(NULL, NULL, "-workdir", wdgrab, sizeof(wdgrab), &flg);
    wd0 = "./";
    if (flg) {
        wd0.append(wdgrab);
        wd0.append("/");
    }


    // Setup sitffness matrix, load vector and bcs (Dirichlet) for the design problem
    SetUpBC(da_nodes,boundaries);

	// Set up optimization retion
	SetUpOptimizationRegion(da_nodes, da_elem);

    // setup timeintergation
    SetUpTimeIntegration(ti);

}

LinearElasticity::~LinearElasticity() {

    // Deallocate
    VecDestroy(&(U0));
    VecDestroy(&(U0_full));
    VecDestroy(&(DU0_full));
    VecDestroy(&(DDU0_full));
    VecDestroy(&(DU0));
    VecDestroy(&(DDU0));
    VecDestroy(&(L0));
    VecDestroy(&(L0_full));
    VecDestroy(&(LT));
    VecDestroy(&(DL0));
    VecDestroy(&(DL0_full));
    VecDestroy(&(DDL0));
    VecDestroy(&(DDL0_full));
    VecDestroy(&(Lp1));
    VecDestroy(&(P));
    VecDestroy(&(DLp1));
    VecDestroy(&(DDLp1));
    VecDestroy(&(R));
    VecDestroyVecs(tsc, &U);
    VecDestroyVecs(tsc, &DU);
    VecDestroyVecs(tsc, &DDU);
    VecDestroyVecs(tsc, &U_full);
    VecDestroyVecs(tsc, &DU_full);
    VecDestroyVecs(tsc, &DDU_full);
    VecDestroyVecs(tsc, &L);
    VecDestroyVecs(tsc, &DL);
    VecDestroyVecs(tsc, &DDL);
    VecDestroyVecs(tsc, &L_full);
    VecDestroyVecs(tsc, &DL_full);
    VecDestroyVecs(tsc, &DDL_full);
    VecDestroy(&RHS);
    VecDestroy(&RHS0);
    VecDestroy(&RHSLocation);
    VecDestroy(&HDU);
    VecDestroy(&HDDU);
    VecDestroy(&HDDDU);
    VecDestroy(&HDU_full);
    VecDestroy(&HDDU_full);
    VecDestroy(&HDDDU_full);
    VecDestroy(&(N));
    VecDestroy(&(Noriginal));
    MatDestroy(&(M_full));
    MatDestroy(&(C_full));
    MatDestroy(&(K_full));
    MatDestroy(&(M));
    MatDestroy(&(C));
    MatDestroy(&(K));
    MatDestroy(&(A));
    MatDestroy(&(dRdS));
    MatDestroy(&(B1));
    MatDestroy(&(B2));
    MatDestroy(&(B3));
    KSPDestroy(&(ksp1));
    KSPDestroy(&(ksp2));
    DMDestroy(&(da_nodal));
	
		if (reductionmodel!=NULL){delete reductionmodel;}

}

PetscErrorCode LinearElasticity::SetUpTimeIntegration(PetscInt ti) {
    PetscErrorCode ierr=0;
    PetscScalar b,g; // Newmark parameters
    switch (ti){
        case 0: // backward euler time integration
        a[0] = 0.0;
        a[1] = 0.0;
        a[2] = 1.0/dt;
        a[3] = 1.0/dt;
        a[4] = 0.0;
        a[5] = 1.0/(dt*dt);
        break;
        case 1:
        b = 0.25;
        g = 0.50;
        a[0] = 1.0-g/b;               // -1
        a[1] = (1.0-g/(2.0*b))*dt;      // zero!
        a[2] =  g/(b*dt);
        a[3] = 1/(b*dt);
        a[4] = 1/(2.0*b)-1.0;           // This!
        a[5] = 1.0/(b*dt*dt);
        break;

    }
    PetscPrintf(PETSC_COMM_WORLD,"#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD,"# T I M E  C O N S T A N T S\n");
    PetscPrintf(PETSC_COMM_WORLD,"# --------------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD,"# a[0]: %f\n",a[0]);
    PetscPrintf(PETSC_COMM_WORLD,"# a[1]: %f\n",a[1]);
    PetscPrintf(PETSC_COMM_WORLD,"# a[2]: %f\n",a[2]);
    PetscPrintf(PETSC_COMM_WORLD,"# a[3]: %f\n",a[3]);
    PetscPrintf(PETSC_COMM_WORLD,"# a[4]: %f\n",a[4]);
    PetscPrintf(PETSC_COMM_WORLD,"# a[5]: %f\n",a[5]);
    return ierr;
}


PetscErrorCode LinearElasticity::SetUpOptimizationRegion(DM da_nodes, DM da_elem){

// INitiaite error code
PetscErrorCode ierr = 0;

// Construct P vector 
ierr = DMCreateGlobalVector(da_elem, &(P));

// Initialize P with ones (default include everything into objective calculation
VecSet(P,0.0);

    PetscScalar     dx, dy, dz;
    DMBoundaryType  bx, by, bz;
    DMDAStencilType stype;

        // Extract information from the nodal mesh
        PetscInt M, N, L, md, nd, pd;
        DMDAGetInfo(da_nodes, NULL, &M, &N, &L, &md, &nd, &pd, NULL, NULL, &bx, &by, &bz, &stype);

        // Find the element size
        Vec lcoor;
        DMGetCoordinatesLocal(da_nodes, &lcoor);
        PetscScalar* lcoorp;
        VecGetArray(lcoor, &lcoorp);

		  // PetscInt        nel, nen;
		  PetscInt nen,nel;
        const PetscInt* necon;
        DMDAGetElements_3D(da_nodes, &nel, &nen, &necon);

        // Use the first element to compute the dx, dy, dz
        dx = lcoorp[3 * necon[0 * nen + 1] + 0] - lcoorp[3 * necon[0 * nen + 0] + 0];
        dy = lcoorp[3 * necon[0 * nen + 2] + 1] - lcoorp[3 * necon[0 * nen + 1] + 1];
        dz = lcoorp[3 * necon[0 * nen + 4] + 2] - lcoorp[3 * necon[0 * nen + 0] + 2];

        nn[0] = M;
        nn[1] = N;
        nn[2] = L;

        ne[0] = nn[0] - 1;
        ne[1] = nn[1] - 1;
        ne[2] = nn[2] - 1;

        xc[0] = 0.0;
        xc[1] = ne[0] * dx;
        xc[2] = 0.0;
        xc[3] = ne[1] * dy;
        xc[4] = 0.0;
        xc[5] = ne[2] * dz;

			
			PetscScalar ex, ey, ez;
			for (PetscInt i = 0; i < nel ; i++){

				  // Find element center coordinate
				  ex = lcoorp[3 * necon[i * nen + 1] + 0]-dx/2;
				  ey = lcoorp[3 * necon[i * nen + 2] + 1]-dy/2;
				  ez = lcoorp[3 * necon[i * nen + 4] + 2]-dz/2;

			
					// USe this to include all elements in optimization 	
			 		//VecSetValueLocal(P, i, 1.00, INSERT_VALUES);	

					switch (optimizationRegion){
					case 0:
					if (ez > xc[5]-1 && ex < xc[0]+1){
					 		VecSetValueLocal(P, i, 1.0, INSERT_VALUES);	
					}
					break;
					case 1:
					if (ez > xc[5]/2.0-0.5 && ez<xc[5]/2.0+0.5 && ex < xc[1]-1){
					 		VecSetValueLocal(P, i, 1.0, INSERT_VALUES);	
					}
					break;
					case 2:
							// include all elements 
					 		VecSetValueLocal(P, i, 1.0, INSERT_VALUES);	
					break;
					case 3:
					//PetscPrintf(PETSC_COMM_WORLD,"element: %i, ex: %f, ey: %f, ez: %f\n", i,ex,ey,ez);
					if (ex > 0.75*xc[1]-2.5 && ex < 0.75*xc[1]+2.5 && 
						 ey > 0.50*xc[3]-2.5 && ey < 0.50*xc[3]+2.5 && 
						 ez > 0.50*xc[5]-2.5 && ez < 0.50*xc[5]+2.5
						 ){
					 		VecSetValueLocal(P, i, 1.0, INSERT_VALUES);	

							PetscPrintf(PETSC_COMM_WORLD,"%f, %f, %f\n",ex,ey,ez);
					}
					break;
					}
				
			}
	
		
	// Compute sum of P
	PetscScalar vecsum;
	VecSum(P,&vecsum);


PetscPrintf(PETSC_COMM_WORLD,"# The objective function is evaluated within %f elements\n",vecsum);
PetscPrintf(PETSC_COMM_WORLD,"# \n",vecsum);


    // Assemble the N vectors
    VecAssemblyBegin(P);
    VecAssemblyEnd(P);

    // restore array 
	 VecRestoreArray(lcoor, &lcoorp);

// return 
return ierr; 

}

PetscErrorCode LinearElasticity::SetUpBC(DM da_nodes,PetscInt boundaries) {

    PetscErrorCode ierr=0;
    // Extract information from input DM and create one for the linear elasticity
    // number of nodal dofs: (u,v,w)
    PetscInt numnodaldof = 3;

    // // Stencil width: each node connects to a box around it - linear elements
    PetscInt stencilwidth = 1;

    PetscScalar     dx, dy, dz;
    DMBoundaryType  bx, by, bz;
    DMDAStencilType stype;
	 PetscInt nel;
    {
        // Extract information from the nodal mesh
        PetscInt M, N, P, md, nd, pd;
        DMDAGetInfo(da_nodes, NULL, &M, &N, &P, &md, &nd, &pd, NULL, NULL, &bx, &by, &bz, &stype);

        // Find the element size
        Vec lcoor;
        DMGetCoordinatesLocal(da_nodes, &lcoor);
        PetscScalar* lcoorp;
        VecGetArray(lcoor, &lcoorp);

		  // PetscInt        nel, nen;
		  PetscInt nen;
        const PetscInt* necon;
        DMDAGetElements_3D(da_nodes, &nel, &nen, &necon);

        // Use the first element to compute the dx, dy, dz
        dx = lcoorp[3 * necon[0 * nen + 1] + 0] - lcoorp[3 * necon[0 * nen + 0] + 0];
        dy = lcoorp[3 * necon[0 * nen + 2] + 1] - lcoorp[3 * necon[0 * nen + 1] + 1];
        dz = lcoorp[3 * necon[0 * nen + 4] + 2] - lcoorp[3 * necon[0 * nen + 0] + 2];

        VecRestoreArray(lcoor, &lcoorp);

        nn[0] = M;
        nn[1] = N;
        nn[2] = P;

        ne[0] = nn[0] - 1;
        ne[1] = nn[1] - 1;
        ne[2] = nn[2] - 1;

        xc[0] = 0.0;
        xc[1] = ne[0] * dx;
        xc[2] = 0.0;
        xc[3] = ne[1] * dy;
        xc[4] = 0.0;
        xc[5] = ne[2] * dz;

    }

    // Create the nodal mesh
    DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, stype, nn[0], nn[1], nn[2], PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, numnodaldof, stencilwidth, 0, 0, 0, &(da_nodal));
    // Initialize
    DMSetFromOptions(da_nodal);
    DMSetUp(da_nodal);

    // Set the coordinates
    DMDASetUniformCoordinates(da_nodal, xc[0], xc[1], xc[2], xc[3], xc[4], xc[5]);
    // Set the element type to Q1: Otherwise calls to GetElements will change to
    // P1 ! STILL DOESN*T WORK !!!!
    DMDASetElementType(da_nodal, DMDA_ELEMENT_Q1);

    // Allocate matrix and the RHS and Solution vector and Dirichlet vector
	 if (reduction ) {
        ierr = DMCreateMatrix(da_nodal, &(M_full));
        MatDuplicate(M_full,MAT_DO_NOT_COPY_VALUES,&C_full);
        MatDuplicate(M_full,MAT_DO_NOT_COPY_VALUES,&K_full);

        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&M);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&C);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&K);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&A);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&B1);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&B2);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&B3);
        MatCreateSeqDense(PETSC_COMM_SELF,basis,basis,NULL,&dRdS);
	 } else
	 {
        ierr = DMCreateMatrix(da_nodal, &(M));
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&C);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&K);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&M_full);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&C_full);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&K_full);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&A);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&B1);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&B2);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&B3);
        MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&dRdS);
	 }
    
	 // Single column vectors
	 if (reduction){
        ierr = DMCreateGlobalVector(da_nodal, &(Noriginal));
        ierr = DMCreateGlobalVector(da_nodal, &(RHSLocation));
        VecCreateSeq(PETSC_COMM_SELF,basis,&N);
        VecCreateSeq(PETSC_COMM_SELF,basis,&U0);
        VecDuplicate(Noriginal,&(U0_full));
        VecDuplicate(Noriginal,&(DU0_full));
        VecDuplicate(Noriginal,&(DDU0_full));
        VecDuplicate(Noriginal,&(L0_full));
        VecDuplicate(Noriginal,&(DL0_full));
        VecDuplicate(Noriginal,&(DDL0_full));
        VecCreateSeq(PETSC_COMM_SELF,basis,&DU0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&DDU0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&L0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&LT);
        VecCreateSeq(PETSC_COMM_SELF,basis,&DL0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&DDL0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&R);
        VecCreateSeq(PETSC_COMM_SELF,basis,&Lp1);
        VecCreateSeq(PETSC_COMM_SELF,basis,&DLp1);
        VecCreateSeq(PETSC_COMM_SELF,basis,&DDLp1);
        VecCreateSeq(PETSC_COMM_SELF,basis,&RHS0);
        VecCreateSeq(PETSC_COMM_SELF,basis,&RHS);
        VecCreateSeq(PETSC_COMM_SELF,basis,&HDU);
        VecCreateSeq(PETSC_COMM_SELF,basis,&HDDU);
        VecCreateSeq(PETSC_COMM_SELF,basis,&HDDDU);
        VecDuplicate(Noriginal,&(HDU_full));
        VecDuplicate(Noriginal,&(HDDU_full));
        VecDuplicate(Noriginal,&(HDDDU_full));
	 } else {
        ierr = DMCreateGlobalVector(da_nodal, &(Noriginal));
        VecDuplicate(Noriginal,&(N));
        VecDuplicate(N,&(U0));
        VecDuplicate(N,&(DU0));
        VecDuplicate(N,&(DDU0));
        VecDuplicate(N,&(L0));
        VecDuplicate(N,&(L0_full));
        VecDuplicate(N,&(U0_full));
        VecDuplicate(N,&(DU0_full));
        VecDuplicate(N,&(DDU0_full));
        VecDuplicate(N,&(LT));
        VecDuplicate(N,&(DL0));
        VecDuplicate(N,&(DDL0));
        VecDuplicate(N,&(DL0_full));
        VecDuplicate(N,&(DDL0_full));
        VecDuplicate(N,&(R));
        VecDuplicate(N,&(Lp1));
        VecDuplicate(N,&(DLp1));
        VecDuplicate(N,&(DDLp1));
        VecDuplicate(N,&(RHS));
        VecDuplicate(N,&(RHS0));
        VecDuplicate(N,&(RHSLocation));
        VecDuplicate(N,&(HDU));
        VecDuplicate(N,&(HDDU));
        VecDuplicate(N,&(HDDDU));
        VecDuplicate(N,&(HDU_full));
        VecDuplicate(N,&(HDDU_full));
        VecDuplicate(N,&(HDDDU_full));
	}
	// MultiColumn vectors
    VecDuplicateVecs(N,tsc, &(U));
    VecDuplicateVecs(N,tsc, &(DU));
    VecDuplicateVecs(N,tsc, &(DDU));
    VecDuplicateVecs(Noriginal,tsc, &(U_full));
    VecDuplicateVecs(Noriginal,tsc, &(DU_full));
    VecDuplicateVecs(Noriginal,tsc, &(DDU_full));
    VecDuplicateVecs(N,tsc, &(L));
    VecDuplicateVecs(N,tsc, &(DL));
    VecDuplicateVecs(N,tsc, &(DDL));
    VecDuplicateVecs(Noriginal,tsc, &(L_full));
    VecDuplicateVecs(Noriginal,tsc, &(DL_full));
    VecDuplicateVecs(Noriginal,tsc, &(DDL_full));

	 if (reduction) {
		 reductionmodel = new SOMMG(Noriginal,basis);
	 }

    // Set the local stiffness matrix
    PetscScalar X[8] = {0.0, dx, dx, 0.0, 0.0, dx, dx, 0.0};
    PetscScalar Y[8] = {0.0, 0.0, dy, dy, 0.0, 0.0, dy, dy};
    PetscScalar Z[8] = {0.0, 0.0, 0.0, 0.0, dz, dz, dz, dz};

    // Compute the element stiffnes matrix - constant due to structured grid
    Hex8IsoparametricKE(X, Y, Z, nu, false, KE);
    Hex8IsoparametricME(X, Y, Z, false, ME);
	PetscPrintf(PETSC_COMM_WORLD,"DALPHA: %f, DBETA: %f\n", dalpha,dbeta);
    Hex8IsoparametricCE(dalpha,dbeta,CE);

    // Set the RHS and Dirichlet vector
    VecSet(U0,0.0);
	 LoadSolutionAsU0();
    VecSet(DU0,0.0);

    // Set the vectors to zero for initial. Filled out later 
    VecSet(RHS, 0.0);
    VecSet(Noriginal, 1.0);
    for (PetscInt j = 0; j<tsc; j++){
        VecSet(U[j], 0.0);
        VecSet(DU[j], 0.0);
        VecSet(DDU[j], 0.0);
    }

    // Global coordinates and a pointer
    Vec          lcoor; // borrowed ref - do not destroy!
    PetscScalar* lcoorp;

    // Get local coordinates in local node numbering including ghosts
    ierr = DMGetCoordinatesLocal(da_nodal, &lcoor);
    CHKERRQ(ierr);
    VecGetArray(lcoor, &lcoorp);

    // Get local dof number
    PetscInt NN;
    VecGetSize(lcoor, &NN);

    // Compute epsilon parameter for finding points in space:
    epsi = PetscMin(dx * 0.05, PetscMin(dy * 0.05, dz * 0.05));

    // Set the values:
    switch (boundaries){
        case 0: // MBB type of boundary conditions 
              for (PetscInt i = 0; i < NN; i++) {
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[0]) < epsi ) {
                      VecSetValueLocal(Noriginal, i, 0.0, INSERT_VALUES);
                      VecSetValueLocal(Noriginal, i+1, 0.0, INSERT_VALUES);
                 }
		           if (i % 3 == 0 && PetscAbsScalar(lcoorp[i]-xc[1]) <epsi && PetscAbsScalar(lcoorp[i+2]-xc[4]) <epsi) {
                      VecSetValueLocal(Noriginal, i+1, 0.0, INSERT_VALUES);
                      VecSetValueLocal(Noriginal, i+2, 0.0, INSERT_VALUES);
		           }
              }
        break;
        case 1: // Cantilever type of boundary conditions
              for (PetscInt i = 0; i < NN; i++) {
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[0]) < epsi ) {
                      VecSetValueLocal(Noriginal, i, 0.0, INSERT_VALUES);
                      VecSetValueLocal(Noriginal, i+1, 0.0, INSERT_VALUES);
                      VecSetValueLocal(Noriginal, i+2, 0.0, INSERT_VALUES);
                 }
              }
        break;
    }


    // Assemble the N vectors
    VecAssemblyBegin(N);
    VecAssemblyEnd(N);
    VecAssemblyBegin(Noriginal);
    VecAssemblyEnd(Noriginal);
	 
    // Boundary conditions on initial conditions
    VecPointwiseMult(U0, U0, N);
    VecPointwiseMult(DU0, DU0, N);

    // Vector restore 
	 VecRestoreArray(lcoor, &lcoorp);

    return ierr;
}

PetscErrorCode LinearElasticity::Load(PetscScalar t) {

    PetscErrorCode ierr=0;

	 // MBB LOAD EXAMPLE MULTI FREQUENCY 
   PetscScalar magnitude = 1;
   PetscScalar t0 = 8; // Time at peak intentisy
   PetscScalar c = 1.1; //sharpness of load. cmaller numbers means narrower load
	PetscScalar omega = 5;//1.5;
	PetscScalar omega1 = 4;
	PetscScalar omega2 = 5;
	PetscScalar omega3 = 6;
	PetscScalar omega4 = 7;
	PetscScalar loadsignal = magnitude*exp(-PetscPowScalar(t-t0,2)/(PetscPowScalar(2*c,2)))*PetscCosReal(omega1*(t-t0))*PetscCosReal(omega2*(t-t0))*PetscCosReal(omega3*(t-t0))*PetscCosReal(omega4*(t-t0));

    // PetscPrintf(PETSC_COMM_WORLD,"t: %f, load is: %f\n",t, loadsignal);
    PetscScalar LoadIntensity = loadsignal/(nn[1]-1);
   
	VecCopy(RHS0,RHS);
	VecScale(RHS,LoadIntensity);

    return ierr;
}
PetscErrorCode LinearElasticity::LoadLocation() {

    PetscErrorCode ierr=0;


    // Global coordinates and a pointer
    Vec          lcoor; // borrowed ref - do not destroy!
    PetscScalar* lcoorp;

    // Get local coordinates in local node numbering including ghosts
    ierr = DMGetCoordinatesLocal(da_nodal, &lcoor);
    CHKERRQ(ierr);
    VecGetArray(lcoor, &lcoorp);

    // Get local dof number
    PetscInt NN;
    VecGetSize(lcoor, &NN);

        //PetscScalar t = dt*j;
        for (PetscInt i = 0; i < NN; i++) {
            switch (loadloc){
				case 0: //MBB-type loading location 
				     // Line load
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[0]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]) < epsi) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0, INSERT_VALUES);
                 }
                 // Adjust the corners
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[0]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]) < epsi && PetscAbsScalar(lcoorp[i + 1] - xc[2]) < epsi ) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0 / 2.0, INSERT_VALUES);
                 }
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[0]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]) < epsi && PetscAbsScalar(lcoorp[i + 1] - xc[3]) < epsi) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0 / 2.0, INSERT_VALUES);
                 }
		       break;
				 case 1: // CANTilever-type load location
				     // Line load
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[1]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]/2.0) < epsi) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0, INSERT_VALUES);
                 }
                 // Adjust the corners
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[1]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]/2.0) < epsi && PetscAbsScalar(lcoorp[i + 1] - xc[2]) < epsi ) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0 / 2.0, INSERT_VALUES);
                 }
                 if (i % 3 == 0 && PetscAbsScalar(lcoorp[i] - xc[1]) < epsi && PetscAbsScalar(lcoorp[i + 2] - xc[5]/2.0) < epsi && PetscAbsScalar(lcoorp[i + 1] - xc[3]) < epsi) {
                     VecSetValueLocal(RHSLocation, i + 2, 1.0 / 2.0, INSERT_VALUES);
                 }
				 break;
				 }

        }
    //}
    VecAssemblyBegin(RHSLocation);
    VecAssemblyEnd(RHSLocation);
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);

    VecRestoreArray(lcoor, &lcoorp);

    return ierr;
}

//PetscErrorCode LinearElasticity::SolveLagrange0(PetscInt c) {
PetscErrorCode LinearElasticity::SolveLagrange0(PetscInt c,Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta) {
    PetscErrorCode ierr=0;

    PetscInt    niter;
    PetscScalar rnorm;
    PetscReal RHSnorm;


    KSPConvergedReason KSPreason;

    //time1 = MPI_Wtime();

    Vec CDDL,KDDL,B1Ltemp,B2Ltemp,B3Ltemp,R1,R2,R3,RR1,RR2,RR3,MM1,MM2,MM3;
    VecDuplicate(R, &CDDL);
    VecDuplicate(R, &KDDL);
    VecDuplicate(R, &B1Ltemp);
    VecDuplicate(R, &B2Ltemp);
    VecDuplicate(R, &B3Ltemp);
    VecDuplicate(R, &R1);
    VecDuplicate(R, &R2);
    VecDuplicate(R, &R3);
    VecDuplicate(R, &RR1);
    VecDuplicate(R, &RR2);
    VecDuplicate(R, &RR3);
    VecDuplicate(R, &MM1);
    VecDuplicate(R, &MM2);
    VecDuplicate(R, &MM3);



    Vec NI;
    VecDuplicate(N, &NI);
	 if (reduction){
	 VecSet(NI,0.0);
	 } else {
    VecSet(NI, 1.0);
    VecAXPY(NI, -1.0, N);
	}


    // reset vectors
    VecSet(R1,0.0);
    VecSet(R2,0.0);
    VecSet(R3,0.0);
    VecSet(B1Ltemp,0.0);
    VecSet(B2Ltemp,0.0);
    

    VecSet(HDU,0.0);
    VecSet(HDDU,0.0);
    VecSet(HDDDU,0.0);
    VecSet(HDU_full,0.0);
    VecSet(HDDU_full,0.0);
    VecSet(HDDDU_full,0.0);

	 // Evaluate gradient of objective wrt the state 
	  ComputeGradientObjectiveElementLoop(xPhys, Emin, Emax,Rmin,Rmax, penalK, penalM,dalpha, dbeta, U0_full, DU0_full,DDU0_full); 

    MatMult(B1,Lp1,B1Ltemp);
    MatMult(B2,Lp1,B2Ltemp);
    MatMult(B3,Lp1,B3Ltemp);

    // Construct residual
    VecCopy(HDU,R1);
    VecAXPY(R1,1.0,B1Ltemp);
    VecAXPY(R1,a[2],DLp1);
    VecAXPY(R1,a[5],DDLp1);
    VecScale(R1,-1.0);

    VecCopy(HDDU,R2);
    VecAXPY(R2,1.0,B2Ltemp);
    VecAXPY(R2,-a[0],DLp1);
    VecAXPY(R2,+a[3],DDLp1);
    VecScale(R2,-1.0);

    VecCopy(HDDDU,R3);
    VecAXPY(R3,1.0,B3Ltemp);
    VecAXPY(R3,-a[1],DLp1);
    VecAXPY(R3,+a[4],DDLp1);
    VecScale(R3,-1.0);

    // BOundary conditions on RHS
    VecPointwiseMult(R3, R3, N);

    // Kopier system matrice til dRdS
    MatCopy(M,dRdS,DIFFERENT_NONZERO_PATTERN);

    // Sæt randbetingelser på dRdS
    MatDiagonalScale(dRdS, N, N);
    MatDiagonalSet(dRdS, NI, ADD_VALUES);

    ierr = KSPSetOperators(ksp, dRdS, dRdS); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);                    CHKERRQ(ierr);

    // Copy last SOLUTION
    VecCopy(DDLp1,DDL0);

	 // Solve for lagrange multipliers 
	 ierr = KSPSolve(ksp, R3, DDL0);

	// GEt DEBUG INFO
    KSPGetIterationNumber(ksp, &niter);
    KSPGetResidualNorm(ksp, &rnorm);
    KSPGetConvergedReason(ksp,&KSPreason);

    VecNorm(R3, NORM_2, &RHSnorm);
    rnorm = rnorm / RHSnorm;

	// Clean up
    VecDestroy(&NI);
    VecDestroy(&R1);
    VecDestroy(&R2);
    VecDestroy(&R3);
    VecDestroy(&RR1);
    VecDestroy(&RR2);
    VecDestroy(&RR3);
    VecDestroy(&MM1);
    VecDestroy(&MM2);
    VecDestroy(&MM3);
    VecDestroy(&CDDL);
    VecDestroy(&KDDL);
    VecDestroy(&B1Ltemp);
    VecDestroy(&B2Ltemp);
    VecDestroy(&B3Ltemp);

    return ierr;
}



PetscErrorCode LinearElasticity::SolveLagrangen(PetscInt c,Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM, PetscScalar dalpha, PetscScalar dbeta) {
    PetscErrorCode ierr=0;

    PetscInt    niter;
    PetscScalar rnorm;
    PetscReal RHSnorm;


    KSPConvergedReason KSPreason;

    //time1 = MPI_Wtime();

    Vec CDDL,KDDL,B1Ltemp,B2Ltemp,B3Ltemp,R1,R2,R3;
    VecDuplicate(R, &CDDL);
    VecDuplicate(R, &KDDL);
    VecDuplicate(R, &B1Ltemp);
    VecDuplicate(R, &B2Ltemp);
    VecDuplicate(R, &B3Ltemp);
    VecDuplicate(R, &R1);
    VecDuplicate(R, &R2);
    VecDuplicate(R, &R3);

    Vec NI;
    VecDuplicate(N, &NI);
	 if (reduction){
	 VecSet(NI,0.0);
	 } else {
    VecSet(NI, 1.0);
    VecAXPY(NI, -1.0, N);
	}

    PetscInt t;

    for (PetscInt j=tsc-1;j>=0;j--){

        // Timestep identifier
        t = tsc*c+j+1;

        // reset vectors
        VecSet(R1,0.0);
        VecSet(R2,0.0);
        VecSet(R3,0.0);
        VecSet(B1Ltemp,0.0);
        VecSet(B2Ltemp,0.0);
        VecSet(B3Ltemp,0.0);

        VecSet(HDU,0.0);
        VecSet(HDDU,0.0);
        VecSet(HDDDU,0.0);
        VecSet(HDU_full,0.0);
        VecSet(HDDU_full,0.0);
        VecSet(HDDDU_full,0.0);

		  ComputeGradientObjectiveElementLoop(xPhys, Emin, Emax,Rmin,Rmax, penalK, penalM,dalpha, dbeta, U_full[j], DU_full[j],DDU_full[j]); 


        if (t == Ts-1 ){

            VecCopy(HDDU,DL[j]);
            VecScale(DL[j],-1.0);

            VecCopy(HDDDU,DDL[j]);
            VecScale(DDL[j],-1.0);

            // Construct residual
            VecCopy(HDU,R1);
            VecAXPY(R1,-a[2],DL[j]);
            VecAXPY(R1,-a[5],DDL[j]);
            VecScale(R1,-1.0);

            // BOundary conditions on RHS
            VecPointwiseMult(R1, R1, N);

            // Retrieve starting guess for backwards solve 
            VecCopy(LT,L[j]);

            // Solve for terminal values of L
            ierr = KSPSolve(ksp, R1, L[j]);

            // Save this solution as the final values of L, for use as starting guess for next backwards solve.
            VecCopy(L[j],LT);

        } else {

            MatMult(B1,Lp1,B1Ltemp);
            MatMult(B2,Lp1,B2Ltemp);
            MatMult(B3,Lp1,B3Ltemp);

            VecCopy(HDDU,DL[j]);
            VecAXPY(DL[j],+1.0,B2Ltemp);
            VecAXPY(DL[j],-a[0],DLp1);
            VecAXPY(DL[j],+a[3],DDLp1);
            VecScale(DL[j],-1.0);

            VecCopy(HDDDU,DDL[j]);
            VecAXPY(DDL[j],1.0,B3Ltemp);
        
            VecAXPY(DDL[j],+a[4],DDLp1);
            VecScale(DDL[j],-1.0);

            // Construct residual
            VecCopy(HDU,R1);
            VecAXPY(R1,+1.0,B1Ltemp);
            VecAXPY(R1,+a[2],DLp1);
            VecAXPY(R1,+a[5],DDLp1);
            VecScale(R1,-1.0);
            VecAXPY(R1,+a[2],DL[j]);
            VecAXPY(R1,+a[5],DDL[j]);

            // Copy last SOLUTION
            VecCopy(Lp1,L[j]);

            // BOundary conditions on RHS
            VecPointwiseMult(R1, R1, N);

				 // SOLVE
            ierr = KSPSolve(ksp, R1, L[j]);

        }

		// Grab debug information 
		KSPGetIterationNumber(ksp, &niter);
      KSPGetResidualNorm(ksp, &rnorm);
      KSPGetConvergedReason(ksp,&KSPreason);
      VecNorm(R1, NORM_2, &RHSnorm);
      rnorm = rnorm / RHSnorm;

	   // Propagate solution	
      VecCopy(L[j],Lp1);
      VecCopy(DL[j],DLp1);
      VecCopy(DDL[j],DDLp1);

    }

	// Clean up
    VecDestroy(&NI);
    VecDestroy(&R1);
    VecDestroy(&R2);
    VecDestroy(&R3);
    VecDestroy(&CDDL);
    VecDestroy(&KDDL);
    VecDestroy(&B1Ltemp);
    VecDestroy(&B2Ltemp);
    VecDestroy(&B3Ltemp);

    return ierr;
}

PetscErrorCode LinearElasticity::CompleteInitialState(Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar Rmin, PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta,PetscBool write,MPIIO *output) {
    PetscErrorCode ierr=0;
    std::string wd = wd0;


    Vec RKU, RCDU, RMDDU;
    VecDuplicate(R, &RKU);
    VecDuplicate(R, &RCDU);
    VecDuplicate(R, &RMDDU);

	 // Setup solver 
	 // Potential multigrid
    if (ksp1 == NULL){
        ierr = SetUpSolver1();
        CHKERRQ(ierr);
    }
	 //For solving reduced systems only LU with cholesky factorization
    if (ksp2 == NULL){
        ierr = SetUpSolver2();
        CHKERRQ(ierr);
    }
	 if (reduction) {
		ksp=ksp2;
	 } else {
		ksp=ksp1;
	 }

    // Assemble the stiffness matrix
    AssembleMassMatrix(xPhys, Rmin, Rmax, penalM);
    AssembleStiffnessMatrix(xPhys, Emin, Emax, penalK);
	 PetscScalar RHSnorm;
    AssembleDampingMatrix(dalpha,dbeta);
    LoadLocation();




	// In case a reduction is enabled
   if (reduction) {
	// Build reduction  basis 
		reductionmodel->ConstructBasis(K_full,M_full, C_full, RHSLocation, Noriginal,omega0, ksp1);
		MatCopy(reductionmodel->Mr,M,SAME_NONZERO_PATTERN);
		MatCopy(reductionmodel->Cr,C,SAME_NONZERO_PATTERN);
		MatCopy(reductionmodel->Kr,K,SAME_NONZERO_PATTERN);
		VecCopy(reductionmodel->Pr,RHS0);
		// Ignore boundary conditions 
		VecSet(N,1.0); 
	} else {
		// Use Noriginal for boundary conditions 
		VecCopy(Noriginal,N);
		VecCopy(RHSLocation,RHS0);
		MatCopy(M_full,M,SAME_NONZERO_PATTERN);
		MatCopy(C_full,C,SAME_NONZERO_PATTERN);
		MatCopy(K_full,K,SAME_NONZERO_PATTERN);
	}

    // Build external force
    Vec NI;
    VecDuplicate(N, &NI);
	 if (reduction){
	 VecSet(NI,0.0);
	 } else {
    VecSet(NI, 1.0);
    VecAXPY(NI, -1.0, N);
	}

    MatZeroEntries(A);
    MatAXPY(A,-a[5],M,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatAXPY(A,-a[2],C,DIFFERENT_NONZERO_PATTERN);// quasi static: 0
    MatAXPY(A,-1.00,K,DIFFERENT_NONZERO_PATTERN);
    MatDiagonalScale(A, N, N);

    // // Create the B1 matrix
    MatZeroEntries(B1);
    MatAXPY(B1,a[5],M,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatAXPY(B1,a[2],C,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatDiagonalScale(B1, N, N);

    // // Create the B2 matrix
    MatZeroEntries(B2);
    MatAXPY(B2,+a[3],M,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatAXPY(B2,-a[0],C,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatDiagonalScale(B2, N, N);

    // // Create the B3 matrix
    MatZeroEntries(B3);
    MatAXPY(B3,+a[4],M,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatAXPY(B3,-a[1],C,DIFFERENT_NONZERO_PATTERN); // quasi static: 0
    MatDiagonalScale(B3, N, N);

    KSPConvergedReason KSPreason;
    PetscInt    niter;
    PetscScalar rnorm;

    // Build
    MatMult(C,DU0,RCDU);
    VecPointwiseMult(RCDU, RCDU, N);

    MatMult(K,U0,RKU);
    VecPointwiseMult(RKU, RKU, N);

    // Build external force
    Load(0.0);

    VecCopy(RHS,R);       // External forcing ??
    VecAXPY(R,-1.0,RCDU);
    VecAXPY(R,-1.0,RKU);

    // Kopier system matrice til dRdS
    MatCopy(M,dRdS,DIFFERENT_NONZERO_PATTERN);

    // Sæt randbetingelser på dRdS
    MatDiagonalScale(dRdS, N, N);
    MatDiagonalSet(dRdS, NI, ADD_VALUES);

    // Randbetingelser på residual
    VecPointwiseMult(R, R, N);

    // set ksp operators and aplpy PC
    ierr = KSPSetOperators(ksp, dRdS, dRdS); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);                    CHKERRQ(ierr);

    ierr = KSPSolve(ksp, R, DDU0);
		
    // get debug information
    KSPGetIterationNumber(ksp, &niter);
    KSPGetResidualNorm(ksp, &rnorm);
    KSPGetConvergedReason(ksp,&KSPreason);
    VecNorm(R, NORM_2, &RHSnorm);
    rnorm = rnorm / RHSnorm;

    // Put initial state to putput if WRITING
    if (write && output != NULL){
		if (reduction) {
			// Part two: Project back to full solution from reduced sol.
			// extract values from bvec
			PetscScalar *pbvec1, *pbvec2, *pbvec3;
			VecGetArray(U0,&pbvec1);
			VecGetArray(DU0,&pbvec2);
			VecGetArray(DDU0,&pbvec3);
			// Zero the fullvec
			VecSet(U0_full, 0.0);
			VecSet(DU0_full, 0.0);
			VecSet(DDU0_full, 0.0);
			for (PetscInt k=0;k<basis;k++){
				VecAXPY(U0_full,pbvec1[k],reductionmodel->QQ[k]);
				VecAXPY(DU0_full,pbvec2[k],reductionmodel->QQ[k]);
				VecAXPY(DDU0_full,pbvec3[k],reductionmodel->QQ[k]);
			}
			VecRestoreArray(U0_full,&pbvec1);
			VecRestoreArray(DU0_full,&pbvec2);
			VecRestoreArray(DDU0_full,&pbvec3);
		}else{
			VecCopy(U0,U0_full);
			VecCopy(DU0,DDU0_full);
			VecCopy(DDU0,DDU0_full);
		}
     	 output->WriteVTK(da_nodal,U0_full,xPhys,0);
			 //}
    }

    // Prepare for battle
    MatCopy(A,dRdS,DIFFERENT_NONZERO_PATTERN);

    // Sæt randbetingelser på dRdS
    MatDiagonalScale(dRdS, N, N);
    MatDiagonalSet(dRdS, NI, ADD_VALUES);

    // complete setup and change of ksp operator, aplpy PC
    ierr = KSPSetOperators(ksp, dRdS, dRdS); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    // Construct filename
    wd = wd.append("cp0000.dat");

    // write output
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,wd.c_str(),FILE_MODE_WRITE,&viewer);
    VecView(U0,viewer);
    VecView(DU0,viewer);
    VecView(DDU0,viewer);
    PetscViewerDestroy(&viewer);


PetscBool dumpToMatlabDebug = PETSC_FALSE;
if (dumpToMatlabDebug) {
PetscPrintf(PETSC_COMM_WORLD,"DUMPING MATRICES FOR DEBUGGIN\n");
        PetscViewer viewer;
		  Mat Mout, Cout,Kout;
		  
    MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&Mout);
    MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&Cout);
    MatDuplicate(M,MAT_DO_NOT_COPY_VALUES,&Kout);

    MatCopy(M,Mout,DIFFERENT_NONZERO_PATTERN);
    MatCopy(C,Cout,DIFFERENT_NONZERO_PATTERN);
    MatCopy(K,Kout,DIFFERENT_NONZERO_PATTERN);

    MatDiagonalScale(Mout, N, N);
    //MatDiagonalSet(Mout, NI, ADD_VALUES);
    MatDiagonalScale(Cout, N, N);
    //MatDiagonalSet(Cout, NI, ADD_VALUES);
    MatDiagonalScale(Kout, N, N);
    MatDiagonalSet(Kout, NI, ADD_VALUES);

		wd = wd0;
	 wd = wd.append("Koutload.m");
		  PetscObjectSetName((PetscObject)Kout,"Kout"); // fortæller at matricen hedder K i matlab
		  PetscViewerASCIIOpen(PETSC_COMM_WORLD, wd.c_str(), &viewer);
		  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
		  MatView(Kout,viewer);
		  PetscViewerPopFormat(viewer);
PetscPrintf(PETSC_COMM_WORLD,"-Done dumping Kout\n");


		wd = wd0;
	 wd = wd.append("Coutload.m");
		  PetscObjectSetName((PetscObject)Cout,"Cout"); // fortæller at matricen hedder K i matlab
		  PetscViewerASCIIOpen(PETSC_COMM_WORLD, wd.c_str(), &viewer);
		  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
		  MatView(Cout,viewer);
		  PetscViewerPopFormat(viewer);
PetscPrintf(PETSC_COMM_WORLD,"-Done dumping Cout\n");

		wd = wd0;
	 wd = wd.append("Moutload.m");
		  PetscObjectSetName((PetscObject)Mout,"Mout"); // fortæller at matricen hedder K i matlab
		  PetscViewerASCIIOpen(PETSC_COMM_WORLD, wd.c_str(), &viewer);
		  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
		  MatView(Mout,viewer);
		  PetscViewerPopFormat(viewer);
PetscPrintf(PETSC_COMM_WORLD,"-Done dumping Mout\n");

PetscViewerDestroy(&viewer);

    MatDestroy(&(Mout));
    MatDestroy(&(Cout));
    MatDestroy(&(Kout));
}


    // time for clearning
    VecDestroy(&NI);
    VecDestroy(&RKU);
    VecDestroy(&RCDU);
    VecDestroy(&RMDDU);

    // Return that errcode
    return ierr;

}

PetscErrorCode LinearElasticity::GenerateCheckPoints(PetscBool write,Vec xPhys,MPIIO * output){
    PetscErrorCode ierr=0;
    char filename[PETSC_MAX_PATH_LEN];
    std::string wd;

    for (PetscInt c = 0; c<ncp; c++){

        // Reset work path
        wd = wd0;

		  PetscReal RHSnorm;
        // Solve the next block
        SolveStateBlock(U0,DU0,DDU0,c);

        if (write && ncp>1 && output != NULL){
            for (PetscInt j = 0; j<tsc; j++){
               // output->WriteVTK(da_nodal,U[j],DU[j],DDU[j],xPhys, tsc*c+j+1);
            if (reduction) {
                // Part two: Project back to full solution from reduced sol.
                // extract values from bvec
                PetscScalar *pbvec1, *pbvec2, *pbvec3;
                VecGetArray(U[j],&pbvec1);
                VecGetArray(DU[j],&pbvec2);
                VecGetArray(DDU[j],&pbvec3);
                // Zero the fullvec
                VecSet(U_full[j], 0.0);
                VecSet(DU_full[j], 0.0);
                VecSet(DDU_full[j], 0.0);
                for (PetscInt k=0;k<basis;k++){
                    VecAXPY(U_full[j],pbvec1[k],reductionmodel->QQ[k]);
                    VecAXPY(DU_full[j],pbvec2[k],reductionmodel->QQ[k]);
                    VecAXPY(DDU_full[j],pbvec3[k],reductionmodel->QQ[k]);
                }
                VecRestoreArray(U_full[j],&pbvec1);
                VecRestoreArray(DU_full[j],&pbvec2);
                VecRestoreArray(DDU_full[j],&pbvec3);
            }else{
                VecCopy(U[j],U_full[j]);
                VecCopy(DU[j],DDU_full[j]);
                VecCopy(DDU[j],DDU_full[j]);
            }
            
                output->WriteVTK(da_nodal,U_full[j],xPhys,tsc*c+j+1);
            
            }
        //			PetscPrintf(PETSC_COMM_WORLD,"Write succeeded.. \n");
        }

        // reset U0 to last of77 U
        VecCopy(U[tsc-1],U0);
        VecCopy(DU[tsc-1],DU0);
        VecCopy(DDU[tsc-1],DDU0);

        // construct filename for checkpoint
        sprintf(filename, "cp%04d.dat", (unsigned short)c+1);

        // create work path with -workdir
        wd = wd.append(filename);

        // check string is correct
        // PetscPrintf(PETSC_COMM_WORLD," - %s\n",wd.c_str());

        // Write checkpoint
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,wd.c_str(),FILE_MODE_WRITE,&viewer);
        VecView(U0,viewer);
        VecView(DU0,viewer);
        VecView(DDU0,viewer);
        PetscViewerDestroy(&viewer);
    }
	 //PetscPrintf(PETSC_COMM_WORLD,"----\n");
    return ierr;
}


PetscErrorCode LinearElasticity::VolumeConstraint(PetscScalar* gx,Vec dgdx, Vec xPhys,PetscScalar volfrac) {

    PetscErrorCode ierr=0;

    // Compute volume constraint gx[0]
        PetscInt neltot;
        VecGetSize(xPhys, &neltot);
        gx[0] = 0;
        VecSum(xPhys, &(gx[0]));

        gx[0] = gx[0] / (((PetscScalar)neltot)) - volfrac;
        VecSet(dgdx, 1.0 / (((PetscScalar)neltot)));
        return ierr;
}

PetscErrorCode LinearElasticity::SolveStateBlock(Vec U0, Vec DU0, Vec DDU0,PetscInt c) {
    PetscErrorCode ierr=0;

    Vec B1Utemp,B2DUtemp,B3DDUtemp,AUtemp,dU,Mddu,Cdu,Ku;
    VecDuplicate(R, &dU);
    VecDuplicate(R, &AUtemp);
    VecDuplicate(R, &B1Utemp);
    VecDuplicate(R, &B2DUtemp);
    VecDuplicate(R, &B3DDUtemp);

    VecDuplicate(R, &Mddu);
    VecDuplicate(R, &Cdu);
    VecDuplicate(R, &Ku);
    
	 KSPConvergedReason KSPreason;
    PetscInt    niter;
    PetscScalar rnorm;
    PetscReal RHSnorm;

    for (PetscInt j = 0; j<tsc; j++){
        // Incrmenet time counter
        //t++;

        t = tsc*c+j+1;

        VecSet(U[j],0.0);
        VecSet(DU[j],0.0);
        VecSet(DDU[j],0.0);

        // Create RHS
        Load(t*dt);

        if (j == 0){
            UM1 = U0;
            DUM1 = DU0;
            DDUM1 =  DDU0;
        }else{
            UM1 = U[j-1];
            DUM1 = DU[j-1];
            DDUM1 = DDU[j-1];
        }


        // Fincd residual contributions
        MatMult(A,U[j],AUtemp); // For linear problems this term is probably always zero?
        MatMult(B1,UM1,B1Utemp); // quasi static: 0
        MatMult(B2,DUM1,B2DUtemp); // quasi static: 0
        MatMult(B3,DDUM1,B3DDUtemp); // quasi static: 0

        // Construct R
        VecCopy(RHS,R);
        VecAXPY(R,+1.0,B1Utemp); // quasi static: 0
        VecAXPY(R,+1.0,B2DUtemp); // quasi static: 0
        VecAXPY(R,+1.0,B3DDUtemp); // quasi static: 0
			VecScale(R,-1.00);

        // SÆT RANDBETINGELSER IGEN
        VecPointwiseMult(R, R, N);

        // Copy previous SOLUTION
        VecCopy(UM1,U[j]);
		  
		  // Solve for displacements
        ierr = KSPSolve(ksp, R, U[j]);

        // Find change in U
        VecCopy(U[j],dU);
        VecAXPY(dU,-1.0,UM1);

        // Find velocities
        VecAXPY(DU[j],+a[0],DUM1);
        VecAXPY(DU[j],+a[1],DDUM1);
        VecAXPY(DU[j],+a[2],dU);

        // Find accelerations
        VecAXPY(DDU[j],-a[3],DUM1);
        VecAXPY(DDU[j],-a[4],DDUM1);
        VecAXPY(DDU[j],+a[5],dU);

        // DEBUG
        // Get iteration number and residual from KSP
        KSPGetIterationNumber(ksp, &niter);
        KSPGetResidualNorm(ksp, &rnorm);
        KSPGetConvergedReason(ksp,&KSPreason);
        VecNorm(R, NORM_2, &RHSnorm);
        rnorm = rnorm / RHSnorm;

    }

    VecDestroy(&dU);
    VecDestroy(&Mddu);
    VecDestroy(&Cdu);
    VecDestroy(&Ku);
    VecDestroy(&AUtemp);
    VecDestroy(&B1Utemp);
    VecDestroy(&B2DUtemp);
    VecDestroy(&B3DDUtemp);

    return ierr;
}


PetscErrorCode LinearElasticity::EvaluateObjectiveScale(PetscScalar xp){

// Errorcode init
PetscErrorCode ierr = 0;

// NO scaling
 s = 1.0; ds = 0.0; 
// Scale objective by xPhys
//  s = xp; ds = 1;
// Scale objective by invrsely with xPhys
// s = 1.0/(xp+1e-9); ds = -1.0/((xp+1e-9)*(xp+1e-9));

return ierr;

}

PetscErrorCode LinearElasticity::EvaluateObjective( PetscScalar kfac, PetscScalar uKu, PetscScalar mfac,PetscScalar duMdu,PetscScalar u,PetscScalar du, PetscScalar uu,PetscScalar dudu,PetscScalar ddu, PetscScalar dduddu,PetscScalar duCdu){

// Errorcode init
PetscErrorCode ierr = 0;

// Compute the objective function [UPDATE THESE]
// TOTAL ENERGY
// h0 = dt*0.5*(kfac * uKu +  mfac * duMdu);
// POTENTIAL ENERGY
//   h0 = dt*0.5*(kfac * uKu);
// KINETIC ENERGY
//  h0 = dt*0.5*(mfac * duMdu);
//SQUARED DISPLACEMENTS
//  h0 = dt*0.5*uu;     // OK
//SQUARED VELOCITIED 
//  h0 = dt*0.5*dudu; // OK 
// SQUARED ACCELERATIONS 
// h0 = dt*0.5*dduddu;
// Return errorcode
//return ierr;



    switch (objective){
        case 0: // squared displacements 
				h0 = dt*0.5*uu;
        break;
        case 1: // squared velocities 
            h0 = dt*0.5*dudu; 
        break;
		  case 2: // squared accelrations
				h0 = dt*0.5*dduddu;
		  break;
		  case 3: // Potential energy
			   h0 = dt*0.5*(kfac * uKu);
		  break;
		  case 4: // Kinetic energy 
				h0 = dt*0.5*(mfac * duMdu);
		  break;
		  case 5: // Total energy 
				h0 = dt*0.5*(kfac * uKu +  mfac * duMdu);
		  break;
		  case 6: // Dissipated energy 
				h0 = dt*0.5*duCdu;
		  break;
    }

	return ierr;

}

// TO HERE

PetscErrorCode LinearElasticity::EvaluateObjectiveGradient( const PetscScalar dkfac, const PetscScalar uKu, const PetscScalar dmfac, const PetscScalar duMdu,const PetscScalar duCdu, const PetscScalar DxduCdu){

// Errorcode init
PetscErrorCode ierr = 0;

// Gradient of the objective wrt physical densities [UPDATE THESE]
// TOTAL ENERGY
// dfx0 = dt*0.5*(dkfac * uKu + dmfac *  duMdu);
// POTENTIAL ENERGY
// dfx0 = dt*0.5*(dkfac * uKu);
// KINETIC ENERGY
//  dfx0 = dt*0.5*(dmfac * duMdu);
// SQUARED DISPLACEMENT
// dfx0 = 0;   // OK
// SQUARED VELOCITIES
//   dfx0 = 0; //OK
// SQUARED ACCELERATIONS 
//    dfx0 = 0; //OK

    switch (objective){
        case 0: // squared displacements 
            dfx0 = 0;
        break;
        case 1: // squared velocities 
            dfx0 = 0;
        break;
		  case 2: // squared accelrations
		      dfx0 = 0;
		  break;
		  case 3: // Potential energy
				dfx0 = dt*0.5*(dkfac * uKu);
		  break;
		  case 4: // Kinetic energy 
				dfx0 = dt*0.5*(dmfac * duMdu);
		  break;
		  case 5: // Total energy 
				dfx0 = dt*0.5*(dkfac * uKu + dmfac *  duMdu);
		  break;
		  case 6: // Dissipated  energy 
				dfx0 = dt*0.5*(DxduCdu);
		  break;
    }
// Return errorcode
return ierr;

}

PetscErrorCode LinearElasticity::GradObjectiveWrtState(PetscInt k, PetscScalar dt,PetscScalar *up,PetscScalar *dup,PetscScalar *ddup,PetscInt *edof,PetscScalar kfac,PetscScalar mfac,PetscScalar dalpha, PetscScalar dbeta){
					 
PetscErrorCode ierr = 0;

// UPDATE THE FOLLOWING. THE OBJECTIVE WRT TO THE STATE

// Reset 
hdu0 = 0.0;
hddu0 = 0.0;
hdddu0 = 0.0;

    switch (objective){
        case 0: // squared displacements 
				hdu0 = dt*up[edof[k]]; 
				hddu0 = 0;
				hdddu0 = 0;
        break;
        case 1: // squared velocities 
				hdu0 = 0; 
				hddu0 = dt*dup[edof[k]];
				hdddu0 = 0;
        break;
		  case 2: // squared accelrations
				hdu0 = 0; 
				hddu0 = 0;
				hdddu0 = dt*ddup[edof[k]];
		  break;
		  case 3: // Potential energy
				for (PetscInt h = 0; h < 24; h++) {
					hdu0+=dt*kfac*KE[k*24+h]*up[edof[h]];
					hddu0+=0;
					hdddu0+=0;
				}
		  break;
		  case 4: // Kinetic energy 
				for (PetscInt h = 0; h < 24; h++) {
					hdu0+=0;
					hddu0+=dt*mfac*ME[k*24+h]*dup[edof[h]];
					hdddu0+=0;
				}
		  break;
		  case 5: // Total energy 
				for (PetscInt h = 0; h < 24; h++) {
					hdu0+=dt*kfac*KE[k*24+h]*up[edof[h]];
					hddu0+=dt*mfac*ME[k*24+h]*dup[edof[h]];
					hdddu0+=0;
				}
		  break;
		  case 6: // dissipated  energy 
				for (PetscInt h = 0; h < 24; h++) {
					hdu0+=0;
					hddu0+=dt*(dalpha*mfac*ME[k*24+h]+dbeta*kfac*KE[k*24+h])*dup[edof[h]];
					hdddu0+=0;
				}
		  break;
    }


return ierr;

}



PetscErrorCode LinearElasticity::InterpolateMass(PetscScalar Rmax, PetscScalar Rmin, PetscScalar penalM, PetscScalar xp){

	// Error code
	PetscErrorCode ierr = 0;

    switch (massInterpolation){
        case 0: //LInear Interpolation
				mfac =  Rmin+xp*(Rmax - Rmin);
				dmfac = Rmax-Rmin;
        break;
        case 1: // SIMP interpolation
				mfac = Rmin + PetscPowScalar(xp, penalM) * (Rmax - Rmin);
				dmfac = penalM * PetscPowScalar(xp, penalM-1) * (Rmax - Rmin);
        break;
		  case 2: // RAMP interpolation
				mfac = Rmin + xp/(1+penalM*(1-xp))*(Rmax-Rmin);
				dmfac = (1+penalM)/PetscPowScalar(1+penalM*(1-xp),2)*(Rmax-Rmin);
		  break;
    }
	 return ierr;
}

PetscErrorCode LinearElasticity::InterpolateStiffness(PetscScalar Emax, PetscScalar Emin, PetscScalar penalK, PetscScalar xp){
	
	// Error code
	PetscErrorCode ierr = 0;

   switch (stiffnessInterpolation){
        case 0: //LInear Interpolation
				kfac =  Emin+xp*(Emax - Emin);
				dkfac = Emax-Emin;
        break;
        case 1: // SIMP interpolation
				kfac = Emin + PetscPowScalar(xp, penalK) * (Emax - Emin);
//				PetscPrintf(PETSC_COMM_WORLD,"# kfac: %f = %f + (%f - %f )*%f^%f\n",kfac,Emin,Emax,Emin,xp[i],penalK);
				dkfac = penalK * PetscPowScalar(xp, penalK-1) * (Emax - Emin);
        break;
		  case 2: // RAMP interpolation
				kfac = Emin + xp/(1+penalK*(1-xp))*(Emax-Emin);
				dkfac = (1+penalK)/PetscPowScalar(1+penalK*(1-xp),2)*(Emax-Emin);
		  break;
   }
	return ierr;
}

PetscErrorCode LinearElasticity::ComputeObjectiveSensitivities(PetscScalar* fx,Vec dfdx,Vec xPhys,  PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, PetscBool write ,MPIIO * output ) {
    // Errorcode
    PetscErrorCode ierr=0;
    std::string wd ;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);

    // Get Solution
    Vec Uloc, DUloc,DDUloc,Ubloc, DUbloc,DDUbloc,Lloc,DLloc,DDLloc;
    DMCreateLocalVector(da_nodal, &Uloc);
    DMCreateLocalVector(da_nodal, &DUloc);
    DMCreateLocalVector(da_nodal, &DDUloc);
    DMCreateLocalVector(da_nodal, &Ubloc);
    DMCreateLocalVector(da_nodal, &DUbloc);
    DMCreateLocalVector(da_nodal, &DDUbloc);
    DMCreateLocalVector(da_nodal, &Lloc);
    DMCreateLocalVector(da_nodal, &DLloc);
    DMCreateLocalVector(da_nodal, &DDLloc);


		// Zero the gradients 
	  VecSet(dfdx,0.0);

// checkpoint filenames
    char filename[PETSC_MAX_PATH_LEN];

    fx[0] = 0.0;
    for (PetscInt c=ncp-1;c>=0;c--) {

        // reset workpath
        wd = wd0;

        // construct filename for checkpoint
        sprintf(filename, "cp%04d.dat", (unsigned short)c);

        // create work path with -workdir
        wd = wd.append(filename);

        //PetscPrintf(PETSC_COMM_WORLD,"---start read, %s\n",filename);
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD,wd.c_str(),FILE_MODE_READ,&viewer);
        VecLoad(U0,viewer);
        VecLoad(DU0,viewer);
        VecLoad(DDU0,viewer);
        PetscViewerDestroy(&viewer);

			// Find correct time
        t = tsc*c;

		// Solve state for this time-block
		SolveStateBlock(U0, DU0, DDU0,c);
		
	// Ensure the full solution is available for otuput and adjoint solve 
	PetscScalar RHSnorm;
    for (PetscInt j = 0; j<tsc; j++){
		if (reduction) {
			// Part two: Project back to full solution from reduced sol.
			// extract values from bvec
			PetscScalar *pbvec1, *pbvec2, *pbvec3;
			VecGetArray(U[j],&pbvec1);
			VecGetArray(DU[j],&pbvec2);
			VecGetArray(DDU[j],&pbvec3);
			// Zero the fullvec
			VecSet(U_full[j], 0.0);
			VecSet(DU_full[j], 0.0);
			VecSet(DDU_full[j], 0.0);
			for (PetscInt k=0;k<basis;k++){
				VecAXPY(U_full[j],pbvec1[k],reductionmodel->QQ[k]);
				VecAXPY(DU_full[j],pbvec2[k],reductionmodel->QQ[k]);
				VecAXPY(DDU_full[j],pbvec3[k],reductionmodel->QQ[k]);
			}
			VecRestoreArray(U_full[j],&pbvec1);
			VecRestoreArray(DU_full[j],&pbvec2);
			VecRestoreArray(DDU_full[j],&pbvec3);
		}else{
		 VecCopy(U[j],U_full[j]);
		 VecCopy(DU[j],DU_full[j]);
		 VecCopy(DDU[j],DDU_full[j]);
		 }
		}


		// Obtain correct U0, DU0 and DDU0 due to checkpointing 
		if (reduction) {
			// Part two: Project back to full solution from reduced sol.
			// extract values from bvec
			PetscScalar *pbvec1, *pbvec2, *pbvec3;
			VecGetArray(U0,&pbvec1);
			VecGetArray(DU0,&pbvec2);
			VecGetArray(DDU0,&pbvec3);
			// Zero the fullvec
			VecSet(U0_full, 0.0);
			VecSet(DU0_full, 0.0);
			VecSet(DDU0_full, 0.0);
			for (PetscInt k=0;k<basis;k++){
				VecAXPY(U0_full,pbvec1[k],reductionmodel->QQ[k]);
				VecAXPY(DU0_full,pbvec2[k],reductionmodel->QQ[k]);
				VecAXPY(DDU0_full,pbvec3[k],reductionmodel->QQ[k]);
			}
			VecRestoreArray(U0_full,&pbvec1);
			VecRestoreArray(DU0_full,&pbvec2);
			VecRestoreArray(DDU0_full,&pbvec3);
		} else {
			VecCopy(U0,U0_full);
			VecCopy(DU0,DU0_full);
			VecCopy(DDU0,DDU0_full);
		}

		// Generate output 
		if (write && ncp==1 && output!=NULL) {
		    //ioutput->WriteVTK(da_nodal,U0_full,xPhys,0);
			   for (PetscInt j = 0; j<tsc; j++){
				   output->WriteVTK(da_nodal,U_full[j],xPhys,tsc*c+j+1);
				}
		}

		// Solve lagrange multipiers 
		ierr = SolveLagrangen(c,xPhys,Emin, Emax, Rmin, Rmax, penalK,  penalM, dalpha, dbeta); 


	//PetscScalar RHSnorm;
    for (PetscInt j = 0; j<tsc; j++){
		if (reduction) {
			// Part two: Project back to full solution from reduced sol.
			// extract values from bvec
			PetscScalar *pbvec1, *pbvec2, *pbvec3;
			VecGetArray(L[j],&pbvec1);
			VecGetArray(DL[j],&pbvec2);
			VecGetArray(DDL[j],&pbvec3);
			// Zero the fullvec
			VecSet(L_full[j], 0.0);
			VecSet(DL_full[j], 0.0);
			VecSet(DDL_full[j], 0.0);
			for (PetscInt k=0;k<basis;k++){
				VecAXPY(L_full[j],pbvec1[k],reductionmodel->QQ[k]);
				VecAXPY(DL_full[j],pbvec2[k],reductionmodel->QQ[k]);
				VecAXPY(DDL_full[j],pbvec3[k],reductionmodel->QQ[k]);
			}
			VecRestoreArray(L_full[j],&pbvec1);
			VecRestoreArray(DL_full[j],&pbvec2);
			VecRestoreArray(DDL_full[j],&pbvec3);
		}else{
		 VecCopy(L[j],L_full[j]);
		 VecCopy(DL[j],DL_full[j]);
		 VecCopy(DDL[j],DDL_full[j]);
		 }
		}


		// COMPUTE OBJECTIVE --------------------------------------------------
        for (PetscInt tc=0; tc<tsc; tc++){
            t = tsc*c+tc+1;
            DMGlobalToLocalBegin(da_nodal, U_full[tc], INSERT_VALUES, Uloc);
            DMGlobalToLocalEnd(da_nodal, U_full[tc], INSERT_VALUES, Uloc);
            DMGlobalToLocalBegin(da_nodal, DU_full[tc], INSERT_VALUES, DUloc);
            DMGlobalToLocalEnd(da_nodal, DU_full[tc], INSERT_VALUES, DUloc);
            DMGlobalToLocalBegin(da_nodal, DDU_full[tc], INSERT_VALUES, DDUloc);
            DMGlobalToLocalEnd(da_nodal, DDU_full[tc], INSERT_VALUES, DDUloc);

				ierr = ComputeObjectiveElementLoop(fx,xPhys, Emin, Emax, Rmin, Rmax, penalK, penalM, dalpha,  dbeta, Uloc, DUloc, DDUloc); 
        }

        // Include initial timestep
        if (c == 0 ){
            t = 0;

            DMGlobalToLocalBegin(da_nodal, U0_full, INSERT_VALUES, Uloc);
            DMGlobalToLocalEnd(da_nodal, U0_full, INSERT_VALUES, Uloc);
            DMGlobalToLocalBegin(da_nodal, DU0_full, INSERT_VALUES, DUloc);
            DMGlobalToLocalEnd(da_nodal, DU0_full, INSERT_VALUES, DUloc);
            DMGlobalToLocalBegin(da_nodal, DDU0_full, INSERT_VALUES, DDUloc);
            DMGlobalToLocalEnd(da_nodal, DDU0_full, INSERT_VALUES, DDUloc);

				ierr = ComputeObjectiveElementLoop(fx, xPhys, Emin, Emax, Rmin, Rmax, penalK, penalM, dalpha,  dbeta, Uloc, DUloc, DDUloc);
        }

		// COMPUTE SENSITIVITIES  -------------------------------------------------
		//
		//
        for (PetscInt tc=0; tc<tsc; tc++){
            t = tsc*c+tc+1;
            DMGlobalToLocalBegin(da_nodal, U_full[tc], INSERT_VALUES, Uloc);
            DMGlobalToLocalEnd(da_nodal, U_full[tc], INSERT_VALUES, Uloc);
            DMGlobalToLocalBegin(da_nodal, DU_full[tc], INSERT_VALUES, DUloc);
            DMGlobalToLocalEnd(da_nodal, DU_full[tc], INSERT_VALUES, DUloc);
            DMGlobalToLocalBegin(da_nodal, DDU_full[tc], INSERT_VALUES, DDUloc);
            DMGlobalToLocalEnd(da_nodal, DDU_full[tc], INSERT_VALUES, DDUloc);

            DMGlobalToLocalBegin(da_nodal, L_full[tc], INSERT_VALUES, Lloc);
            DMGlobalToLocalEnd(da_nodal, L_full[tc], INSERT_VALUES, Lloc);
            DMGlobalToLocalBegin(da_nodal, DL_full[tc], INSERT_VALUES, DLloc);
            DMGlobalToLocalEnd(da_nodal, DL_full[tc], INSERT_VALUES, DLloc);
            DMGlobalToLocalBegin(da_nodal, DDL_full[tc], INSERT_VALUES, DDLloc);
            DMGlobalToLocalEnd(da_nodal, DDL_full[tc], INSERT_VALUES, DDLloc);

            if (tc==0) {
                DMGlobalToLocalBegin(da_nodal, U0_full, INSERT_VALUES, Ubloc);
                DMGlobalToLocalEnd(da_nodal, U0_full, INSERT_VALUES, Ubloc);
                DMGlobalToLocalBegin(da_nodal, DU0_full, INSERT_VALUES, DUbloc);
                DMGlobalToLocalEnd(da_nodal, DU0_full, INSERT_VALUES, DUbloc);
                DMGlobalToLocalBegin(da_nodal, DDU0_full, INSERT_VALUES, DDUbloc);
                DMGlobalToLocalEnd(da_nodal, DDU0_full, INSERT_VALUES, DDUbloc);
            }else{
                DMGlobalToLocalBegin(da_nodal, U_full[tc-1], INSERT_VALUES, Ubloc);
                DMGlobalToLocalEnd(da_nodal, U_full[tc-1], INSERT_VALUES, Ubloc);
                DMGlobalToLocalBegin(da_nodal, DU_full[tc-1], INSERT_VALUES, DUbloc);
                DMGlobalToLocalEnd(da_nodal, DU_full[tc-1], INSERT_VALUES, DUbloc);
                DMGlobalToLocalBegin(da_nodal, DDU_full[tc-1], INSERT_VALUES, DDUbloc);
                DMGlobalToLocalEnd(da_nodal, DDU_full[tc-1], INSERT_VALUES, DDUbloc);
            }

	 ierr = ComputeSensitivitiesElementLoopn(dfdx,  xPhys,  Emin, Emax, Rmin, Rmax,  penalK,  penalM, dalpha, dbeta, Uloc, DUloc, DDUloc,Lloc,  Ubloc, DUbloc, DDUbloc,N,P);

        }

        // Include initial timestep
        if (c == 0 ){
            t = 0;
//	PetscPrintf(PETSC_COMM_WORLD,"======================================\n");

            // Solve initial timesteps lagranges
			    ierr = SolveLagrange0(c,xPhys,Emin, Emax, Rmin, Rmax, penalK,  penalM, dalpha, dbeta); 
				if (reduction) {
					// extract values from bvec
					PetscScalar *pbvec1, *pbvec2, *pbvec3;
					VecGetArray(L0,&pbvec1);
					VecGetArray(DL0,&pbvec2);
					VecGetArray(DDL0,&pbvec3);
					// Zero the fullvec
					VecSet(L0_full, 0.0);
					VecSet(DL0_full, 0.0);
					VecSet(DDL0_full, 0.0);
					for (PetscInt k=0;k<basis;k++){
						VecAXPY(L0_full,pbvec1[k],reductionmodel->QQ[k]);
						VecAXPY(DL0_full,pbvec2[k],reductionmodel->QQ[k]);
						VecAXPY(DDL0_full,pbvec3[k],reductionmodel->QQ[k]);
					}
					VecRestoreArray(L0_full,&pbvec1);
					VecRestoreArray(DL0_full,&pbvec2);
					VecRestoreArray(DDL0_full,&pbvec3);
				} 
				else 
				{
					VecCopy(L0,L0_full);
					VecCopy(DL0,DL0_full);
					VecCopy(DDL0,DDL0_full);
				}	

            DMGlobalToLocalBegin(da_nodal, U0_full, INSERT_VALUES, Uloc);
            DMGlobalToLocalEnd(da_nodal, U0_full, INSERT_VALUES, Uloc);
            DMGlobalToLocalBegin(da_nodal, DU0_full, INSERT_VALUES, DUloc);
            DMGlobalToLocalEnd(da_nodal, DU0_full, INSERT_VALUES, DUloc);
            DMGlobalToLocalBegin(da_nodal, DDU0_full, INSERT_VALUES, DDUloc);
            DMGlobalToLocalEnd(da_nodal, DDU0_full, INSERT_VALUES, DDUloc);

            DMGlobalToLocalBegin(da_nodal, DDL0_full, INSERT_VALUES, DDLloc);
            DMGlobalToLocalEnd(da_nodal, DDL0_full, INSERT_VALUES, DDLloc);

            ierr = ComputeSensitivitiesElementLoop0(dfdx,  xPhys,  Emin, Emax, Rmin, Rmax,  penalK,  penalM, dalpha, dbeta, Uloc, DUloc, DDUloc, DDLloc, N,P);
        }

    }



            PetscScalar tmp = fx[0];
            fx[0]           = 0.0;
            MPI_Allreduce(&tmp, &(fx[0]), 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);

            VecDestroy(&Uloc);
            VecDestroy(&DUloc);
            VecDestroy(&DDUloc);
            VecDestroy(&Ubloc);
            VecDestroy(&DUbloc);
            VecDestroy(&DDUbloc);
            VecDestroy(&Lloc);
            VecDestroy(&DLloc);
            VecDestroy(&DDLloc);

    return ierr;
}


PetscErrorCode LinearElasticity::ComputeSensitivitiesElementLoopn( Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc,Vec DDUloc,Vec Lloc,Vec  Ubloc,Vec DUbloc, Vec DDUbloc,Vec Ni,Vec P) {

    // Errorcode
    PetscErrorCode ierr=0;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);

    // Edof array
    PetscInt edof[24];
    //PetscScalar dAdx[546];

    // Get dfdx
    PetscScalar* df;
    PetscScalar dfx;
    VecGetArray(dfdx, &df);

    // Get pointer to the densities
    PetscScalar* xp;
    VecGetArray(xPhys, &xp);

    // get pointer to local vector
    PetscScalar* up;
    PetscScalar* dup;
    PetscScalar* ddup;
    PetscScalar* ubp;
    PetscScalar* dubp;
    PetscScalar* ddubp;
    PetscScalar* lp;
	 PetscScalar* np;
	 PetscScalar* pp;

    VecGetArray(Ubloc, &ubp);
    VecGetArray(DUbloc, &dubp);
    VecGetArray(DDUbloc, &ddubp);
    VecGetArray(Uloc, &up);
    VecGetArray(DUloc, &dup);
    VecGetArray(DDUloc, &ddup);
    VecGetArray(Lloc, &lp);
	 VecGetArray(N,&np);
	 VecGetArray(P,&pp);

    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }
		

			// interpolate mass properties 
			InterpolateMass(Rmax,Rmin,penalM,xp[i]);
			InterpolateStiffness(Emax,Emin,penalK,xp[i]);

        //Terms used for the gradient calculations
        PetscScalar t2 = 0.0;
        PetscScalar t3 = 0.0;
        PetscScalar uKu = 0.0;
        PetscScalar duMdu = 0.0;
        PetscScalar duCdu = 0.0;
		  PetscScalar DxduCdu = 0.0;
		  PetscScalar u = 0.0;
		  PetscScalar du = 0.0;
		  PetscScalar uu = 0.0;
		  PetscScalar dudu = 0.0;
		  PetscScalar ddu = 0.0;
		  PetscScalar dduddu = 0.0;
			PetscScalar dAdxudBdxu = 0.0;
			PetscScalar sdAdxu = 0.0;
			PetscScalar sdBdxu = 0.0;
			PetscScalar st3_1 = 0.0;
			PetscScalar st3_2 = 0.0;
			PetscScalar st3_3 = 0.0;

        for (PetscInt k = 0; k < 24; k++) {
		  //if (i==1) {
			//		PetscPrintf(PETSC_COMM_WORLD,"lp[edof[%i]]: %e\n",k,lp[edof[k]]);
			//		}
            for (PetscInt h = 0; h < 24; h++) {

                PetscScalar dme = dmfac*ME[k * 24 + h]; // OK
                PetscScalar dke = dkfac*KE[k * 24 + h]; // OK
                PetscScalar dce = dalpha*dme + dbeta*dke; // OK

                PetscScalar dAdx  = -dke - a[5]*dme - a[2]*dce; // OK
                PetscScalar dBdx1 = a[5]*dme+a[2]*dce; // OK
                PetscScalar dBdx2 = a[3]*dme-a[0]*dce; // OK
                PetscScalar dBdx3 = a[4]*dme-a[1]*dce; // OK

                t2 = 0.0;
                
					 PetscScalar t3_1 = dBdx1*ubp[edof[h]];
                PetscScalar t3_2 = dBdx2*dubp[edof[h]];
                PetscScalar t3_3 = dBdx3*ddubp[edof[h]];
                PetscScalar dAdxu = dAdx*up[edof[h]];
                PetscScalar dBdxu = t3_1 + t3_2 + t3_3;
               
					//if (i==1){
					//	PetscPrintf(PETSC_COMM_WORLD,"%t3: %4.10f\n",t3);
					//}
					 t3 += lp[edof[k]] * (dAdxu + dBdxu);

					dAdxudBdxu += dAdxu+dBdxu;
					sdAdxu += dAdxu;
					sdBdxu += dBdxu;
					
					st3_1 += t3_1;
					st3_2 += t3_2;
					st3_3 += t3_3;
                
					 uKu += up[edof[k]] * KE[k * 24 + h] * up[edof[h]];
                duMdu += dup[edof[k]] * ME[k * 24 + h] * dup[edof[h]];
                duCdu += dup[edof[k]] * CE[k*24 +h]* dup[edof[h]];
                DxduCdu += dup[edof[k]] * dce * dup[edof[h]];
            }

            // if (i==0){
            //    PetscPrintf(PETSC_COMM_WORLD," t3: %e,  lp[edof[%i]]: %e dBdx1: %e, dBdx2: %e, dBdx3: %e\n",t3,k, lp[edof[k]]);
            // }

				u += up[edof[k]];
				uu += up[edof[k]]*up[edof[k]];
				du += dup[edof[k]];
				dudu += dup[edof[k]]*dup[edof[k]];
				ddu += ddup[edof[k]];
				dduddu += ddup[edof[k]]*ddup[edof[k]];

				//if (i==1){
				//	PetscPrintf(PETSC_COMM_WORLD,"u[%i]=%e\n",k,up[edof[k]]);
				//}
        }
		  //if (i==1){
		  //PetscPrintf(PETSC_COMM_WORLD,"uu = %e\n",uu);
		  //}
       
        // Add to objective
        EvaluateObjective(kfac, uKu, mfac, duMdu,u,du,uu,dudu,ddu,dduddu,duCdu);
		  EvaluateObjectiveGradient(dkfac, uKu, dmfac, duMdu,duCdu,DxduCdu);
		  EvaluateObjectiveScale(xp[i]);

			// Allow for indicators and elementwise scalings 
			dfx = pp[i]*(h0*ds+dfx0*s);
        
			// Set the Senstivity
			df[i]+=dfx+t2+t3;
		
    }


    VecRestoreArray(dfdx, &df);
    VecRestoreArray(xPhys, &xp);
    VecRestoreArray(Uloc, &up);
    VecRestoreArray(DUloc, &dup);
    VecRestoreArray(DDUloc, &ddup);
    VecRestoreArray(Ubloc, &up);
    VecRestoreArray(DUbloc, &dup);
    VecRestoreArray(DDUbloc, &ddup);
    VecRestoreArray(N, &np);
    VecRestoreArray(P, &pp);

    return ierr;
}


PetscErrorCode LinearElasticity::ComputeSensitivitiesElementLoop0( Vec dfdx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc,Vec DDUloc, Vec DDLloc,Vec N, Vec P) {

    // Errorcode
    PetscErrorCode ierr=0;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);

    // Edof array
    PetscInt edof[24];
    //PetscScalar dAdx[546];

    // Get dfdx
    PetscScalar* df;
	 PetscScalar dfx;
    VecGetArray(dfdx, &df);

    // Get pointer to the densities
    PetscScalar* xp;
    VecGetArray(xPhys, &xp);

    // get pointer to local vector
    PetscScalar* up;
    PetscScalar* dup;
    PetscScalar* ddup;
    PetscScalar* ddlp;
	 PetscScalar* np;
	 PetscScalar* pp;

	// Get local arrays 
    VecGetArray(Uloc, &up);
    VecGetArray(DUloc, &dup);
    VecGetArray(DDUloc, &ddup);
    VecGetArray(DDLloc, &ddlp);
    VecGetArray(N, &np);
    VecGetArray(P, &pp);

    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }

			// Interpolate element properties 
			InterpolateMass(Rmax,Rmin,penalM,xp[i]);
			InterpolateStiffness(Emax,Emin,penalK,xp[i]);

        //Terms used for the gradient calculations
        PetscScalar t2 = 0.0;
        PetscScalar t3 = 0.0;
        PetscScalar uKu = 0.0;
        PetscScalar duMdu = 0.0;
		  PetscScalar duCdu = 0.0;
		  PetscScalar DxduCdu = 0.0;
		  PetscScalar u = 0.0;
		  PetscScalar du = 0.0;
		  PetscScalar uu = 0.0;
		  PetscScalar dudu = 0.0;
		  PetscScalar ddu = 0.0;
		  PetscScalar dduddu = 0.0;
			PetscScalar dAdxudBdxu = 0.0;
			PetscScalar sdAdxu = 0.0;
			PetscScalar sdBdxu = 0.0;

        for (PetscInt k = 0; k < 24; k++) {
            for (PetscInt h = 0; h < 24; h++) {

					 // Interpolated element matrices 
                PetscScalar dme = dmfac*ME[k * 24 + h]; // OK
                PetscScalar dke = dkfac*KE[k * 24 + h]; // OK
                PetscScalar dce = dalpha*dme + dbeta*dke; // OK

                // this product should only be evaluated for the inital timestep
                PetscScalar t2_1 = ddlp[edof[k]] * dke * up[edof[h]];
                PetscScalar t2_2 = ddlp[edof[k]] * dce * dup[edof[h]];
                PetscScalar t2_3 = ddlp[edof[k]] * dme * ddup[edof[h]];
                t2 +=  t2_1 +  t2_2 + t2_3;
                 
					// T3 in this timestep  is zero 
					t3 = 0.0;

                uKu += up[edof[k]] * KE[k * 24 + h] * up[edof[h]];
                duMdu += dup[edof[k]] * ME[k * 24 + h] * dup[edof[h]];
                duCdu += dup[edof[k]] * CE[k * 24 + h] * dup[edof[h]];
                DxduCdu += dup[edof[k]] * dce * dup[edof[h]];
            }

				u += up[edof[k]];
				uu += up[edof[k]]*up[edof[k]];
				du += dup[edof[k]];
				dudu += dup[edof[k]]*dup[edof[k]];
				ddu += ddup[edof[k]];
				dduddu += ddup[edof[k]]*ddup[edof[k]];
				
        }
        // Add to objective
        EvaluateObjective( kfac, uKu, mfac, duMdu,u,du,uu,dudu,ddu,dduddu,duCdu);
		  EvaluateObjectiveGradient( dkfac, uKu, dmfac, duMdu,duCdu,DxduCdu);
		  EvaluateObjectiveScale(xp[i]);
				
			// Allow for elementwise scaling 
			dfx = pp[i]*(h0*ds+dfx0*s);

        // Set the Senstivity
        df[i]+=dfx+t2+t3;
    }

	 // Clean up
    VecRestoreArray(dfdx, &df);
    VecRestoreArray(xPhys, &xp);
    VecRestoreArray(Uloc, &up);
    VecRestoreArray(DUloc, &dup);
    VecRestoreArray(DDUloc, &ddup);
    VecRestoreArray(DDLloc, &ddlp);
    VecRestoreArray(N, &np);
    VecRestoreArray(P, &np);

    return ierr;
}

PetscErrorCode LinearElasticity::ComputeGradientObjectiveElementLoop( Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec U, Vec DU,Vec DDU) {


Vec Uloc, DUloc, DDUloc, HDUloc, HDDUloc, HDDDUloc;

    // Errorcode
    PetscErrorCode ierr=0;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);

	DMCreateLocalVector(da_nodal,&Uloc);
	DMCreateLocalVector(da_nodal,&DUloc);
	DMCreateLocalVector(da_nodal,&DDUloc);
	DMCreateLocalVector(da_nodal,&HDUloc);
	DMCreateLocalVector(da_nodal,&HDDUloc);
	DMCreateLocalVector(da_nodal,&HDDDUloc);
	 

	// zero local vectors 
	VecSet(HDUloc,0.0);
	VecSet(HDDUloc,0.0);
	VecSet(HDDDUloc,0.0);


	DMGlobalToLocalBegin(da_nodal,U,INSERT_VALUES,Uloc);
	DMGlobalToLocalEnd(da_nodal,U,INSERT_VALUES,Uloc);
	DMGlobalToLocalBegin(da_nodal,DU,INSERT_VALUES,DUloc);
	DMGlobalToLocalEnd(da_nodal,DU,INSERT_VALUES,DUloc);
	DMGlobalToLocalBegin(da_nodal,DDU,INSERT_VALUES,DDUloc);
	DMGlobalToLocalEnd(da_nodal,DDU,INSERT_VALUES,DDUloc);
	 
    // Edof array
    PetscInt edof[24];
	 PetscScalar hdu[24];
	 PetscScalar hddu[24];
	 PetscScalar hdddu[24];

    // get pointer to local vector
    PetscScalar* up;
    PetscScalar* dup;
    PetscScalar* ddup;
	 PetscScalar* np;
	 PetscScalar* pp;
	 PetscScalar* xp;

    VecGetArray(Uloc, &up);
    VecGetArray(DUloc, &dup);
    VecGetArray(DDUloc, &ddup);
	 VecGetArray(N,&np);
	 VecGetArray(P,&pp);
	 VecGetArray(xPhys,&xp);



    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }

			// Compute objective scaling 
			EvaluateObjectiveScale(xp[i]);
		  
		  // interpolate mass properties 
		  InterpolateMass(Rmax,Rmin,penalM,xp[i]);
		  InterpolateStiffness(Emax,Emin,penalK,xp[i]);



        for (PetscInt k = 0; k < 24; k++) {

			   GradObjectiveWrtState(k,dt,up,dup,ddup,edof,kfac,mfac,dalpha,dbeta);

				// Construct gradient terms 
				hdu[k] = pp[i]*s*(hdu0);
				hddu[k] = pp[i]*s*(hddu0);
				hdddu[k] = pp[i]*s*(hdddu0);

        }

		 // Assemble gradients of objective wrt state 	
		 VecSetValues(HDUloc,24,edof,hdu,ADD_VALUES);
		 VecSetValues(HDDUloc,24,edof,hddu,ADD_VALUES);
		 VecSetValues(HDDDUloc,24,edof,hdddu,ADD_VALUES);
		 

	}



// PetscPrintf(PETSC_COMM_WORLD,"---\n");
 VecAssemblyBegin(HDU_full);
 VecAssemblyEnd(HDU_full);
 VecAssemblyBegin(HDDU_full);
 VecAssemblyEnd(HDDU_full);
 VecAssemblyBegin(HDDDU_full);
 VecAssemblyEnd(HDDDU_full);


// Construct global vectores 
DMLocalToGlobal(da_nodal, HDUloc,ADD_VALUES,HDU_full);
DMLocalToGlobal(da_nodal, HDDUloc,ADD_VALUES,HDDU_full);
DMLocalToGlobal(da_nodal, HDDDUloc,ADD_VALUES,HDDDU_full);


if (reduction) {
			// REDUCE HDU, HDDU, HDDDU
			PetscScalar dotval;
			Vec tmp1,tmp2,tmp3;
			VecDuplicate(HDU,&tmp1);
			VecDuplicate(HDDU,&tmp2);
			VecDuplicate(HDDDU,&tmp3);
			for (PetscInt i=0;i<basis;i++){
				VecDot(HDU_full,reductionmodel->QQ[i],&dotval);
				VecSetValues(HDU,1,&i,&dotval,INSERT_VALUES);
				VecDot(HDDU_full,reductionmodel->QQ[i],&dotval);
				VecSetValues(HDDU,1,&i,&dotval,INSERT_VALUES);
				VecDot(HDDDU_full,reductionmodel->QQ[i],&dotval);
				VecSetValues(HDDDU,1,&i,&dotval,INSERT_VALUES);
			 }
			 VecDestroy(&tmp1);
			 VecDestroy(&tmp2);
			 VecDestroy(&tmp3);
} else {
VecCopy(HDU_full,HDU);
VecCopy(HDDU_full,HDDU);
VecCopy(HDDDU_full,HDDDU);
}


    VecRestoreArray(U, &up);
    VecRestoreArray(DU, &dup);
    VecRestoreArray(DDU, &ddup);
    VecRestoreArray(N, &np);
    VecRestoreArray(P, &pp);
    VecRestoreArray(xPhys, &xp);
VecDestroy(&HDUloc);
VecDestroy(&HDDUloc);
VecDestroy(&HDDDUloc);
VecDestroy(&Uloc);
VecDestroy(&DUloc);
VecDestroy(&DDUloc);
    return ierr;
}


PetscErrorCode LinearElasticity::ComputeObjectiveElementLoop(PetscScalar* fx, Vec xPhys, PetscScalar Emin,PetscScalar Emax,PetscScalar Rmin,PetscScalar Rmax, PetscScalar penalK, PetscScalar penalM,PetscScalar dalpha, PetscScalar dbeta, Vec Uloc, Vec DUloc,Vec DDUloc) {

    // Errorcode
    PetscErrorCode ierr=0;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);

    // Edof array
    PetscInt edof[24];

    // get pointer to local vector
    PetscScalar* up;
    PetscScalar* dup;
    PetscScalar* ddup;
	 PetscScalar* np;
	 PetscScalar* xp;
	 PetscScalar* pp;

    VecGetArray(Uloc, &up);
    VecGetArray(DUloc, &dup);
    VecGetArray(DDUloc, &ddup);
	 VecGetArray(N,&np);
	 VecGetArray(P,&pp);
	 VecGetArray(xPhys,&xp);

    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }

			// Compute objective scaling 
			EvaluateObjectiveScale(xp[i]);

		  // interpolate mass properties 
		  InterpolateMass(Rmax,Rmin,penalM,xp[i]);
		  InterpolateStiffness(Emax,Emin,penalK,xp[i]);

        //Terms used for the gradient calculations
        PetscScalar uKu = 0.0;
        PetscScalar duMdu = 0.0;
		  PetscScalar duCdu = 0.0;
		  PetscScalar DxduCdu = 0.0;
		  PetscScalar u = 0.0;
		  PetscScalar du = 0.0;
		  PetscScalar uu = 0.0;
		  PetscScalar dudu = 0.0;
		  PetscScalar ddu = 0.0;
		  PetscScalar dduddu = 0.0;

        for (PetscInt k = 0; k < 24; k++) {
            
				for (PetscInt h = 0; h < 24; h++) {
					 
					 uKu += up[edof[k]] * KE[k * 24 + h] * up[edof[h]];
                duMdu += dup[edof[k]] * ME[k * 24 + h] * dup[edof[h]];

                duCdu += dup[edof[k]] * (dalpha*mfac*ME[k * 24 + h] + dbeta*kfac*KE[k*24 + h ])* dup[edof[h]];
                DxduCdu += dup[edof[k]] * (dalpha*dmfac*ME[k * 24 + h] + dbeta*dkfac*KE[k*24 + h ])* dup[edof[h]];
            }
			//PetscPrintf(PETSC_COMM_SELF,"EvaluateObjective up: %e\n",up[edof[k]]);		
				u += up[edof[k]];
				uu += up[edof[k]]*up[edof[k]];
				du += dup[edof[k]];
				dudu += dup[edof[k]]*dup[edof[k]];
				ddu += ddup[edof[k]];
				dduddu += ddup[edof[k]]*ddup[edof[k]];

        }
       
        //  Elementwise contributionto objective
        EvaluateObjective(kfac, uKu, mfac, duMdu,u,du,uu,dudu,ddu,dduddu,duCdu);

		  // Sum the objective  
		  fx[0] += h0*s*pp[i];

    }

    VecRestoreArray(Uloc, &up);
    VecRestoreArray(DUloc, &dup);
    VecRestoreArray(DDUloc, &ddup);
    VecRestoreArray(N, &np);
    VecRestoreArray(P, &np);
    VecRestoreArray(xPhys, &xp);

    return ierr;
}




PetscErrorCode LinearElasticity::WriteRestartFiles() {

    PetscErrorCode ierr = 0;

    // Only dump data if correct allocater has been used
    if (!restart) {
        return -1;
    }

    // Choose previous set of restart files
    if (flip) {
        flip = PETSC_FALSE;
    } else {
        flip = PETSC_TRUE;
    }

    // Open viewers for writing
    PetscViewer view; // vectors
    if (!flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename00.c_str(), FILE_MODE_WRITE, &view);
    } else if (flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename01.c_str(), FILE_MODE_WRITE, &view);
    }

	// Vectores needed for a perfect restart
	VecView(DDU0,view);
   VecView(LT,view);
	VecView(DDL0,view);

    // Clean up
    PetscViewerDestroy(&view);
    
	 if (!flip) {
		  PetscPrintf(PETSC_COMM_WORLD,"# Restart file written to 00-files.\n");
    } else if (flip) {
		  PetscPrintf(PETSC_COMM_WORLD,"# Restart file written to 01-files.\n");
    }

    return ierr;
}

//##################################################################
//######################## PRIVATE #################################
//##################################################################

PetscErrorCode LinearElasticity::AssembleStiffnessMatrix(Vec xPhys, PetscScalar min, PetscScalar max,PetscScalar penal) {

    PetscErrorCode ierr=0;

    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);

    // Get pointer to the densities
    PetscScalar* xp;
    VecGetArray(xPhys, &xp);

    // Zero the matrix
    MatZeroEntries(K_full);

    // Edof array
    PetscInt    edof[24];
    PetscScalar ke[24 * 24];

    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }
        // Use SIMP for stiffness interpolation
		  InterpolateStiffness(max,min,penal,xp[i]);

        for (PetscInt k = 0; k < 24 * 24; k++) {
				ke[k] = KE[k] * kfac;

        }
        // Add values to the sparse matrix
        ierr = MatSetValuesLocal(K_full, 24, edof, 24, edof, ke, ADD_VALUES);
        CHKERRQ(ierr);
    }
    MatAssemblyBegin(K_full, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K_full, MAT_FINAL_ASSEMBLY);

    // VecDestroy(&NI);
    VecRestoreArray(xPhys, &xp);
    DMDARestoreElements(da_nodal, &nel, &nen, &necon);

    return ierr;
}
PetscErrorCode LinearElasticity::AssembleMassMatrix(Vec xPhys, PetscScalar min, PetscScalar max,PetscScalar penal) {


    PetscErrorCode ierr=0;
    // Get the FE mesh structure (from the nodal mesh)
    PetscInt        nel, nen;
    const PetscInt* necon;
    ierr = DMDAGetElements_3D(da_nodal, &nel, &nen, &necon);
    CHKERRQ(ierr);

    // Get pointer to the densities
    PetscScalar* xp;
    VecGetArray(xPhys, &xp);

    // Zero the matrix
    MatZeroEntries(M_full);

    // Edof array
    PetscInt    edof[24];
    PetscScalar me[24 * 24];
		//PetscScalar mesum;

    // Loop over elements
    for (PetscInt i = 0; i < nel; i++) {
        // loop over element nodes
        for (PetscInt j = 0; j < nen; j++) {
            // Get local dofs
            for (PetscInt k = 0; k < 3; k++) {
                edof[j * 3 + k] = 3 * necon[i * nen + j] + k;
            }
        }
        // Use SIMP for stiffness interpolation
		 InterpolateMass(max,min,penal,xp[i]);
		 //mesum = 0.0;
        for (PetscInt k = 0; k < 24 * 24; k++) {
			me[k] = ME[k]*mfac;
        }
		  //PetscPrintf(PETSC_COMM_WORLD,"mesum: %e\n",mesum);
        ierr = MatSetValuesLocal(M_full, 24, edof, 24, edof, me, ADD_VALUES);
    }
    MatAssemblyBegin(M_full, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M_full, MAT_FINAL_ASSEMBLY);


    // VecDestroy(&NI);
    VecRestoreArray(xPhys, &xp);
    DMDARestoreElements(da_nodal, &nel, &nen, &necon);



    return ierr; 
}

PetscErrorCode LinearElasticity::AssembleDampingMatrix(PetscScalar dalpha, PetscScalar dbeta) {



    PetscErrorCode ierr=0;

    // Zero the matrix
    MatZeroEntries(C_full);


    ierr = MatAXPY(C_full,dalpha,M_full,SAME_NONZERO_PATTERN);
    ierr = MatAXPY(C_full,dbeta,K_full,SAME_NONZERO_PATTERN);

    return ierr;
}

PetscErrorCode LinearElasticity::LoadSolutionAsU0() {

    PetscErrorCode ierr=0;

    // CHECK FOR RESTART POINT
    PetscBool flg, onlyDesign;
    onlyDesign = PETSC_FALSE;
	 PetscBool solutionAsU0 = PETSC_FALSE;
    char filenameChar[PETSC_MAX_PATH_LEN];
    PetscOptionsGetBool(NULL, NULL, "-onlyLoadDesign", &onlyDesign,&flg); // DONT READ DESIGN IF THIS IS TRUE
    PetscOptionsGetBool(NULL, NULL, "-solutionAsU0", &solutionAsU0,&flg); // DONT READ DESIGN IF THIS IS TRUE

        // THE FILES FOR WRITING RESTARTS
        std::string filenameWorkdir = "./";
        PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);
        if (flg) {
            filenameWorkdir = "";
            filenameWorkdir.append(filenameChar);
        }
        filename00 = filenameWorkdir;
        filename01 = filenameWorkdir;
        filename00.append("/RestartSol00.dat");
        filename01.append("/RestartSol01.dat");

            
				// Where to read the restart point from
            std::string restartFileVec = ""; // NO RESTART FILE !!!!!
            // GET FILENAME
            PetscOptionsGetString(NULL, NULL, "-restartFileVecSol", filenameChar, sizeof(filenameChar), &flg);
            if (flg) {
                restartFileVec.append(filenameChar);
            }

				if (onlyDesign && solutionAsU0) {

            // PRINT TO SCREEN
            PetscPrintf(PETSC_COMM_WORLD,"# Loaded solution from %s as initial displacements U0 (-restartFileVecSol):  \n", restartFileVec.c_str());
                // Check if files exist:
                PetscBool vecFile = fexists(restartFileVec);
                if (!vecFile) {
                    PetscPrintf(PETSC_COMM_WORLD, "# File: %s NOT FOUND \n", restartFileVec.c_str());
                }

                // READ
                if (vecFile) {
                    PetscViewer view;
                    // Open the data files
                    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, restartFileVec.c_str(), FILE_MODE_READ, &view);

								// This vector is available from restart files in the original static code 
								VecLoad(U0,view);
						  
						  PetscViewerDestroy(&view);
                }
            }
	return ierr;				
}

PetscErrorCode LinearElasticity::SetUpSolver1() {

    PetscErrorCode ierr=0;

    // CHECK FOR RESTART POINT
    restart = PETSC_TRUE;
    flip    = PETSC_TRUE;
    PetscBool flg, onlyDesign;
    onlyDesign = PETSC_FALSE;
    char filenameChar[PETSC_MAX_PATH_LEN];
    PetscOptionsGetBool(NULL, NULL, "-restart", &restart, &flg);
    PetscOptionsGetBool(NULL, NULL, "-onlyLoadDesign", &onlyDesign,&flg); // DONT READ DESIGN IF THIS IS TRUE

    // READ THE RESTART FILE INTO THE SOLUTION VECTOR(S)
    if (restart) {
        // THE FILES FOR WRITING RESTARTS
        std::string filenameWorkdir = "./";
        PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);
        if (flg) {
            filenameWorkdir = "";
            filenameWorkdir.append(filenameChar);
        }
        filename00 = filenameWorkdir;
        filename01 = filenameWorkdir;
        filename00.append("/RestartSol00.dat");
        filename01.append("/RestartSol01.dat");

        // CHECK FOR SOLUTION AND READ TO STATE VECTOR(s)
            
				// Where to read the restart point from
            std::string restartFileVec = ""; // NO RESTART FILE !!!!!
            // GET FILENAME
            PetscOptionsGetString(NULL, NULL, "-restartFileVecSol", filenameChar, sizeof(filenameChar), &flg);
            if (flg) {
                restartFileVec.append(filenameChar);
            }

				if (!onlyDesign) {

            // PRINT TO SCREEN
            PetscPrintf(PETSC_COMM_WORLD,"# Restarting with solution (State Vector) from (-restartFileVecSol): %s \n", restartFileVec.c_str());
                // Check if files exist:
                PetscBool vecFile = fexists(restartFileVec);
                if (!vecFile) {
                    PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileVec.c_str());
                }

                // READ
                if (vecFile) {
                    PetscViewer view;
                    // Open the data files
                    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, restartFileVec.c_str(), FILE_MODE_READ, &view);

								VecLoad(DDU0,view);
								VecLoad(LT,view);
								VecLoad(DDL0,view);
                    
						  PetscViewerDestroy(&view);
                }
            }
				
        }

        PC pc;

        // The fine grid Krylov method

        KSPCreate(PETSC_COMM_WORLD, &(ksp1));

        // SET THE DEFAULT SOLVER PARAMETERS
        // The fine grid solver settings
        PetscScalar rtol         = 1.0e-10; //5 // Use e-10 for FD check 
        PetscScalar atol         = 1.0e-50;
        PetscScalar dtol         = 1.0e5;
        PetscInt    restart      = 100;
        PetscInt    maxitsGlobal = 200;

        // Coarsegrid solver
        PetscScalar coarse_rtol    = 1.0e-10; //8 //Use e-10 for FDcheck
        PetscScalar coarse_atol    = 1.0e-50;
        PetscScalar coarse_dtol    = 1e5;
        PetscInt    coarse_maxits  = 30;
        PetscInt    coarse_restart = 30;

        // Number of smoothening iterations per up/down smooth_sweeps
        PetscInt smooth_sweeps = 4;

        // Set up the solver
        ierr = KSPSetType(ksp1, KSPFGMRES); // KSPCG, KSPGMRES
        CHKERRQ(ierr);

        ierr = KSPGMRESSetRestart(ksp1, restart);
        CHKERRQ(ierr);

        ierr = KSPSetTolerances(ksp1, rtol, atol, dtol, maxitsGlobal);
        CHKERRQ(ierr);


	 PetscBool initialGuessNonzero;
    PetscOptionsGetBool(NULL, NULL, "-ksp_initial_guess_nonzero", &initialGuessNonzero, &flg);
			
		  if (initialGuessNonzero) {
        ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_TRUE);
		  } else {
        ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_FALSE);
		  }

        CHKERRQ(ierr);

        // The preconditinoer
        KSPGetPC(ksp1, &pc);
        // Make PCMG the default solver
        PCSetType(pc, PCMG);

        // Set solver from options
        KSPSetFromOptions(ksp1);

        // Get the prec again - check if it has changed
        KSPGetPC(ksp1, &pc);

        // Flag for pcmg pc
        PetscBool pcmg_flag = PETSC_TRUE;
        PetscObjectTypeCompare((PetscObject)pc, PCMG, &pcmg_flag);

        // Only if PCMG is used
        if (pcmg_flag) {

            // DMs for grid hierachy
            DM *da_list, *daclist;
            Mat R;

            PetscMalloc(sizeof(DM) * nlvls, &da_list);
            for (PetscInt k = 0; k < nlvls; k++)
            da_list[k] = NULL;
            PetscMalloc(sizeof(DM) * nlvls, &daclist);
            for (PetscInt k = 0; k < nlvls; k++)
            daclist[k] = NULL;

            // Set 0 to the finest level
            daclist[0] = da_nodal;

            // Coordinates
            PetscReal xmin = xc[0], xmax = xc[1], ymin = xc[2], ymax = xc[3], zmin = xc[4], zmax = xc[5];

            // Set up the coarse meshes
            DMCoarsenHierarchy(da_nodal, nlvls - 1, &daclist[1]);
            for (PetscInt k = 0; k < nlvls; k++) {
                // NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ???
                da_list[k] = daclist[nlvls - 1 - k];
                // THIS SHOULD NOT BE NECESSARY
                DMDASetUniformCoordinates(da_list[k], xmin, xmax, ymin, ymax, zmin, zmax);
            }


            // the PCMG specific options
            PCMGSetLevels(pc, nlvls, NULL);
            PCMGSetType(pc, PC_MG_MULTIPLICATIVE); // Default
            ierr = PCMGSetCycleType(pc, PC_MG_CYCLE_V);
            CHKERRQ(ierr);
            PCMGSetGalerkin(pc, PC_MG_GALERKIN_BOTH);

            for (PetscInt k = 1; k < nlvls; k++) {
                DMCreateInterpolation(da_list[k - 1], da_list[k], &R, NULL);
                PCMGSetInterpolation(pc, k, R);
                MatDestroy(&R);
            }


            // tidy up
            for (PetscInt k = 1; k < nlvls; k++) { // DO NOT DESTROY LEVEL 0
                DMDestroy(&daclist[k]);
            }
            PetscFree(da_list);
            PetscFree(daclist);


            // AVOID THE DEFAULT FOR THE MG PART
            {
                // SET the coarse grid solver:
                // i.e. get a pointer to the ksp and change its settings
                KSP cksp;
                PCMGGetCoarseSolve(pc, &cksp);
                // The solver
                ierr = KSPSetType(cksp, KSPGMRES); // KSPCG, KSPFGMRES
                ierr = KSPGMRESSetRestart(cksp, coarse_restart);

                ierr = KSPSetTolerances(cksp, coarse_rtol, coarse_atol, coarse_dtol, coarse_maxits);
                // The preconditioner
                PC cpc;
                KSPGetPC(cksp, &cpc);
                PCSetType(cpc, PCSOR); // PCGAMG, PCSOR, PCSPAI (NEEDS TO BE COMPILED), PCJACOBI


					// Coarse solver settings  
					ierr = KSPSetFromOptions(cksp);
					CHKERRQ(ierr); 
					ierr = PCSetFromOptions(cpc);
					CHKERRQ(ierr); 

                // Set smoothers on all levels (except for coarse grid):
                for (PetscInt k = 1; k < nlvls; k++) {

                    KSP dksp;
                    PCMGGetSmoother(pc, k, &dksp);
                    PC dpc;
                    KSPGetPC(dksp, &dpc);
                    ierr = KSPSetType(dksp, KSPGMRES); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)

                    ierr = KSPGMRESSetRestart(dksp, smooth_sweeps);
                    // ierr = KSPSetType(dksp,KSPCHEBYSHEV);
                    ierr = KSPSetTolerances(dksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                        smooth_sweeps); // NOTE in the above maxitr=restart;
                        PCSetType(dpc, PCSOR);                  // PCJACOBI, PCSOR for KSPCHEBYSHEV very good
					
					// Coarse solver settings  
					ierr = KSPSetFromOptions(dksp);
					CHKERRQ(ierr); 
					ierr = PCSetFromOptions(dpc);
					CHKERRQ(ierr); 
               
						}

                }
            }


            // Write check to screen:
            // Check the overall Krylov solver
            KSPType ksptype;
            KSPGetType(ksp1, &ksptype);
            PCType pctype;
            PCGetType(pc, &pctype);
            PetscInt mmax;
            KSPGetTolerances(ksp1, NULL, NULL, NULL, &mmax);
            PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
            PetscPrintf(PETSC_COMM_WORLD, "# L I N E A R  S O L V E R  S E T T I N G S\n");
				PetscPrintf(PETSC_COMM_WORLD, "# ------------------------------------------------------------\n");
            PetscPrintf(PETSC_COMM_WORLD, "# Main solver: %s, prec.: %s, maxiter.: %i \n", ksptype, pctype, mmax);

            // Only if pcmg is used
            if (pcmg_flag) {
                // Check the smoothers and coarse grid solver:
                for (PetscInt k = 0; k < nlvls; k++) {
                    KSP     dksp;
                    PC      dpc;
                    KSPType dksptype;
                    PCMGGetSmoother(pc, k, &dksp);
                    KSPGetType(dksp, &dksptype);
                    KSPGetPC(dksp, &dpc);
                    PCType dpctype;
                    PCGetType(dpc, &dpctype);
                    PetscInt mmax;
                    KSPGetTolerances(dksp, NULL, NULL, NULL, &mmax);
                    PetscPrintf(PETSC_COMM_WORLD, "# Level %i smoother: %s, prec.: %s, sweep: %i \n", k, dksptype, dpctype,
                    mmax);
                }
            }
            PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");

            return ierr;
        }
	
PetscErrorCode LinearElasticity::SetUpSolver2() {

    PetscErrorCode ierr=0;

    // CHECK FOR RESTART POINT
    restart = PETSC_TRUE;
    flip    = PETSC_TRUE;
    PetscBool flg, onlyDesign;
    onlyDesign = PETSC_FALSE;
    char filenameChar[PETSC_MAX_PATH_LEN];
    PetscOptionsGetBool(NULL, NULL, "-restart", &restart, &flg);
    PetscOptionsGetBool(NULL, NULL, "-onlyLoadDesign", &onlyDesign,&flg); // DONT READ DESIGN IF THIS IS TRUE

    // READ THE RESTART FILE INTO THE SOLUTION VECTOR(S)
    if (restart) {
        // THE FILES FOR WRITING RESTARTS
        std::string filenameWorkdir = "./";
        PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);
        if (flg) {
            filenameWorkdir = "";
            filenameWorkdir.append(filenameChar);
        }
        filename00 = filenameWorkdir;
        filename01 = filenameWorkdir;
        filename00.append("/RestartSol00.dat");
        filename01.append("/RestartSol01.dat");

        // CHECK FOR SOLUTION AND READ TO STATE VECTOR(s)
            
				// Where to read the restart point from
            std::string restartFileVec = ""; // NO RESTART FILE !!!!!
            // GET FILENAME
            PetscOptionsGetString(NULL, NULL, "-restartFileVecSol", filenameChar, sizeof(filenameChar), &flg);
            if (flg) {
                restartFileVec.append(filenameChar);
            }

				if (!onlyDesign) {

            // PRINT TO SCREEN
            PetscPrintf(PETSC_COMM_WORLD,"# Restarting with solution (State Vector) from (-restartFileVecSol): %s \n", restartFileVec.c_str());
                // Check if files exist:
                PetscBool vecFile = fexists(restartFileVec);
                if (!vecFile) {
                    PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileVec.c_str());
                }

                // READ
                if (vecFile) {
                    PetscViewer view;
                    // Open the data files
                    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, restartFileVec.c_str(), FILE_MODE_READ, &view);

								VecLoad(DDU0,view);
								VecLoad(LT,view);
								VecLoad(DDL0,view);
                    
						  PetscViewerDestroy(&view);
                }
            }
				
        }


        // The fine grid Krylov method

        KSPCreate(PETSC_COMM_SELF, &(ksp2));

        // SET THE DEFAULT SOLVER PARAMETERS
		KSPSetType(ksp2, KSPPREONLY);
		PC pcloc;
		KSPGetPC(ksp2, &pcloc);
		PCSetFromOptions(pcloc); // Note - all runtime settings overwritten next!
		PCSetType(pcloc, PCLU);

            // Write check to screen:
            // Check the overall Krylov solver
            KSPType ksptype;
            KSPGetType(ksp2, &ksptype);
            PCType pctype;
            PCGetType(pcloc, &pctype);
            PetscInt mmax;
            KSPGetTolerances(ksp2, NULL, NULL, NULL, &mmax);
            PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
            PetscPrintf(PETSC_COMM_WORLD, "# L I N E A R  S O L V E R  S E T T I N G S \n");
				PetscPrintf(PETSC_COMM_WORLD, "# ------------------------------------------------------------\n");
            PetscPrintf(PETSC_COMM_WORLD, "# Main solver: %s, prec.: %s, maxiter.: %i \n", ksptype, pctype, mmax);
            PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
				//}
            return ierr;
        }

        PetscErrorCode LinearElasticity::DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]) {
            PetscErrorCode ierr=0;
            DM_DA*         da = (DM_DA*)dm->data;
            PetscInt       i, xs, xe, Xs, Xe;
            PetscInt       j, ys, ye, Ys, Ye;
            PetscInt       k, zs, ze, Zs, Ze;
            PetscInt       cnt = 0, cell[8], ns = 1, nn = 8;
            PetscInt       c;
            if (!da->e) {
                if (da->elementtype == DMDA_ELEMENT_Q1) {
                    ns = 1;
                    nn = 8;
                }
                ierr = DMDAGetCorners(dm, &xs, &ys, &zs, &xe, &ye, &ze);
                CHKERRQ(ierr);
                ierr = DMDAGetGhostCorners(dm, &Xs, &Ys, &Zs, &Xe, &Ye, &Ze);
                CHKERRQ(ierr);
                xe += xs;
                Xe += Xs;
                if (xs != Xs)
                xs -= 1;
                ye += ys;
                Ye += Ys;
                if (ys != Ys)
                ys -= 1;
                ze += zs;
                Ze += Zs;
                if (zs != Zs)
                zs -= 1;
                da->ne = ns * (xe - xs - 1) * (ye - ys - 1) * (ze - zs - 1);
                PetscMalloc((1 + nn * da->ne) * sizeof(PetscInt), &da->e);
                for (k = zs; k < ze - 1; k++) {
                    for (j = ys; j < ye - 1; j++) {
                        for (i = xs; i < xe - 1; i++) {
                            cell[0] = (i - Xs) + (j - Ys) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                            cell[1] = (i - Xs + 1) + (j - Ys) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                            cell[2] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                            cell[3] = (i - Xs) + (j - Ys + 1) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                            cell[4] = (i - Xs) + (j - Ys) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                            cell[5] = (i - Xs + 1) + (j - Ys) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                            cell[6] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                            cell[7] = (i - Xs) + (j - Ys + 1) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                            if (da->elementtype == DMDA_ELEMENT_Q1) {
                                for (c = 0; c < ns * nn; c++)
                                da->e[cnt++] = cell[c];
                            }
                        }
                    }
                }
            }
            *nel = da->ne;
            *nen = nn;
            *e   = da->e;
            return (0);
        }

        PetscInt LinearElasticity::Hex8IsoparametricKE(PetscScalar* X, PetscScalar* Y, PetscScalar* Z, PetscScalar nu,PetscInt redInt, PetscScalar* ke) {
            // HEX8_ISOPARAMETRIC - Computes HEX8 isoparametric element matrices
            // The element stiffness matrix is computed as:
            //
            //       ke = int(int(int(B^T*C*B,x),y),z)
            //
            // For an isoparameteric element this integral becomes:
            //
            //       ke = int(int(int(B^T*C*B*det(J),xi=-1..1),eta=-1..1),zeta=-1..1)
            //
            // where B is the more complicated expression:
            // B = [dx*alpha1 + dy*alpha2 + dz*alpha3]*N
            // where
            // dx = [invJ11 invJ12 invJ13]*[dxi deta dzeta]
            // dy = [invJ21 invJ22 invJ23]*[dxi deta dzeta]
            // dy = [invJ31 invJ32 invJ33]*[dxi deta dzeta]
            //
            // Remark: The elasticity modulus is left out in the below
            // computations, because we multiply with it afterwards (the aim is
            // topology optimization).
            // Furthermore, this is not the most efficient code, but it is readable.
            //
            /////////////////////////////////////////////////////////////////////////////////
            //////// INPUT:
            // X, Y, Z  = Vectors containing the coordinates of the eight nodes
            //               (x1,y1,z1,x2,y2,z2,...,x8,y8,z8). Where node 1 is in the
            //               lower left corner, and node 2 is the next node
            //               counterclockwise (looking in the negative z-dir). Finish the
            //               x-y-plane and then move in the positive z-dir.
            // redInt   = Reduced integration option boolean (here an integer).
            //           	redInt == 0 (false): Full integration
            //           	redInt == 1 (true): Reduced integration
            // nu 		= Poisson's ratio.
            //
            //////// OUTPUT:
            // ke  = Element stiffness matrix. Needs to be multiplied with elasticity
            // modulus
            //
            //   Written 2013 at
            //   Department of Mechanical Engineering
            //   Technical University of Denmark (DTU).
            /////////////////////////////////////////////////////////////////////////////////

            //// COMPUTE ELEMENT STIFFNESS MATRIX
            // Lame's parameters (with E=1.0):
            PetscScalar lambda = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
            PetscScalar mu     = 1.0 / (2.0 * (1.0 + nu));

//				PetscPrintf(PETSC_COMM_WORLD,"lambda: %f\n",lambda);
//				PetscPrintf(PETSC_COMM_WORLD,"mu: %f\n",mu);


            // Constitutive matrix
            PetscScalar C[6][6] = {{lambda + 2.0 * mu, lambda, lambda, 0.0, 0.0, 0.0},
            {lambda, lambda + 2.0 * mu, lambda, 0.0, 0.0, 0.0},
            {lambda, lambda, lambda + 2.0 * mu, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, mu, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, mu, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, mu}};
            // Gauss points (GP) and weigths
            // Two Gauss points in all directions (total of eight)
            PetscScalar GP[2] = {-0.577350269189626, 0.577350269189626};
            // Corresponding weights
            PetscScalar W[2] = {1.0, 1.0};
            // If reduced integration only use one GP
            if (redInt) {
                GP[0] = 0.0;
                W[0]  = 2.0;
            }
            // Matrices that help when we gather the strain-displacement matrix:
            PetscScalar alpha1[6][3];
            PetscScalar alpha2[6][3];
            PetscScalar alpha3[6][3];
            memset(alpha1, 0, sizeof(alpha1[0][0]) * 6 * 3); // zero out
            memset(alpha2, 0, sizeof(alpha2[0][0]) * 6 * 3); // zero out
            memset(alpha3, 0, sizeof(alpha3[0][0]) * 6 * 3); // zero out
            alpha1[0][0] = 1.0;
            alpha1[3][1] = 1.0;
            alpha1[5][2] = 1.0;
            alpha2[1][1] = 1.0;
            alpha2[3][0] = 1.0;
            alpha2[4][2] = 1.0;
            alpha3[2][2] = 1.0;
            alpha3[4][1] = 1.0;
            alpha3[5][0] = 1.0;
            PetscScalar  dNdxi[8];
            PetscScalar  dNdeta[8];
            PetscScalar  dNdzeta[8];
            PetscScalar  J[3][3];
            PetscScalar  invJ[3][3];
            PetscScalar  beta[6][3];
            PetscScalar  B[6][24]; // Note: Small enough to be allocated on stack
            PetscScalar* dN;
            // Make sure the stiffness matrix is zeroed out:
            memset(ke, 0, sizeof(ke[0]) * 24 * 24);
            // Perform the numerical integration
            for (PetscInt ii = 0; ii < 2 - redInt; ii++) {
                for (PetscInt jj = 0; jj < 2 - redInt; jj++) {
                    for (PetscInt kk = 0; kk < 2 - redInt; kk++) {
                        // Integration point
                        PetscScalar xi   = GP[ii];
                        PetscScalar eta  = GP[jj];
                        PetscScalar zeta = GP[kk];
                        // Differentiated shape functions
                        DifferentiatedShapeFunctions(xi, eta, zeta, dNdxi, dNdeta, dNdzeta);
                        // Jacobian
                        J[0][0] = Dot(dNdxi, X, 8);
                        J[0][1] = Dot(dNdxi, Y, 8);
                        J[0][2] = Dot(dNdxi, Z, 8);
                        J[1][0] = Dot(dNdeta, X, 8);
                        J[1][1] = Dot(dNdeta, Y, 8);
                        J[1][2] = Dot(dNdeta, Z, 8);
                        J[2][0] = Dot(dNdzeta, X, 8);
                        J[2][1] = Dot(dNdzeta, Y, 8);
                        J[2][2] = Dot(dNdzeta, Z, 8);
                        // Inverse and determinant
                        PetscScalar detJ = Inverse3M(J, invJ);
                        // Weight factor at this point
                        PetscScalar weight = W[ii] * W[jj] * W[kk] * detJ;
                        // Strain-displacement matrix
                        memset(B, 0, sizeof(B[0][0]) * 6 * 24); // zero out
                        for (PetscInt ll = 0; ll < 3; ll++) {
                            // Add contributions from the different derivatives
                            if (ll == 0) {
                                dN = dNdxi;
                            }
                            if (ll == 1) {
                                dN = dNdeta;
                            }
                            if (ll == 2) {
                                dN = dNdzeta;
                            }
                            // Assemble strain operator
                            for (PetscInt i = 0; i < 6; i++) {
                                for (PetscInt j = 0; j < 3; j++) {
                                    beta[i][j] =
                                    invJ[0][ll] * alpha1[i][j] + invJ[1][ll] * alpha2[i][j] + invJ[2][ll] * alpha3[i][j];
                                }
                            }
                            // Add contributions to strain-displacement matrix
                            for (PetscInt i = 0; i < 6; i++) {
                                for (PetscInt j = 0; j < 24; j++) {
                                    B[i][j] = B[i][j] + beta[i][j % 3] * dN[j / 3];
                                }
                            }
                        }
                        // Finally, add to the element matrix
                        for (PetscInt i = 0; i < 24; i++) {
                            for (PetscInt j = 0; j < 24; j++) {
                                for (PetscInt k = 0; k < 6; k++) {
                                    for (PetscInt l = 0; l < 6; l++) {
                                        ke[j + 24 * i] = ke[j + 24 * i] + weight * (B[k][i] * C[k][l] * B[l][j]);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return 0;
        }

        PetscInt LinearElasticity::Hex8IsoparametricME(PetscScalar* X, PetscScalar* Y, PetscScalar* Z, PetscInt redInt, PetscScalar* me) {
            // HEX8_ISOPARAMETRIC - Computes HEX8 isoparametric element matrices
            // The element mass  matrix is computed as:
            //
            /////////////////////////////////////////////////////////////////////////////////
            //////// INPUT:
            // X, Y, Z  = Vectors containing the coordinates of the eight nodes
            //               (x1,y1,z1,x2,y2,z2,...,x8,y8,z8). Where node 1 is in the
            //               lower left corner, and node 2 is the next node
            //               counterclockwise (looking in the negative z-dir). Finish the
            //               x-y-plane and then move in the positive z-dir.
            // redInt   = Reduced integration option boolean (here an integer).
            //           	redInt == 0 (false): Full integration
            //           	redInt == 1 (true): Reduced integration
            // nu 		= Poisson's ratio.
            //
            //////// OUTPUT:
            // me  = Element mass matrix. Needs to be multiplied with mass-density
            // modulus
            //
            /////////////////////////////////////////////////////////////////////////////////
            // Two Gauss points in all directions (total of eight)
            PetscScalar GP[2] = {-0.577350269189626, 0.577350269189626};
            // Corresponding weights
            PetscScalar W[2] = {1.0, 1.0};
            // If reduced integration only use one GP
            if (redInt) {
                GP[0] = 0.0;
                W[0]  = 2.0;
            }
            PetscScalar  S[8];
            PetscScalar  vS[72];
            PetscScalar  dNdxi[8];
            PetscScalar  dNdeta[8];
            PetscScalar  dNdzeta[8];
            PetscScalar  J[3][3];
            PetscScalar  invJ[3][3];
            //PetscScalar  B[6][24]; // Note: Small enough to be allocated on stack
            // Make sure the mass matrix is zeroed out:
            memset(me, 0, sizeof(me[0]) * 24 * 24);
            // Perform the numerical integration
            for (PetscInt ii = 0; ii < 2 - redInt; ii++) {
                for (PetscInt jj = 0; jj < 2 - redInt; jj++) {
                    for (PetscInt kk = 0; kk < 2 - redInt; kk++) {
                        // Integration point
                        PetscScalar xi   = GP[ii];
                        PetscScalar eta  = GP[jj];
                        PetscScalar zeta = GP[kk];
                        // Shape functions
                        ShapeFunctions(xi, eta, zeta, S);
                        ShapeFunctionsVector(S,vS);

                        // Differentiated shape functions
                        DifferentiatedShapeFunctions(xi, eta, zeta, dNdxi, dNdeta, dNdzeta);
                        // Jacobian
                        J[0][0] = Dot(dNdxi, X, 8);
                        J[0][1] = Dot(dNdxi, Y, 8);
                        J[0][2] = Dot(dNdxi, Z, 8);
                        J[1][0] = Dot(dNdeta, X, 8);
                        J[1][1] = Dot(dNdeta, Y, 8);
                        J[1][2] = Dot(dNdeta, Z, 8);
                        J[2][0] = Dot(dNdzeta, X, 8);
                        J[2][1] = Dot(dNdzeta, Y, 8);
                        J[2][2] = Dot(dNdzeta, Z, 8);
                        // Inverse and determinant
                        PetscScalar detJ = Inverse3M(J, invJ);
                        // Weight factor at this point
                        PetscScalar weight = W[ii] * W[jj] * W[kk] * detJ;
                        // Finally, add to the element matrix
                        for (PetscInt i = 0; i < 24; i++) {
                            for (PetscInt j = 0; j < 24; j++) {
                                for (PetscInt k = 0; k < 3; k++){
                                    me[j + 24 * i] = me[j + 24 * i] + weight*(vS[i+k*24] * vS[j+24*k]);
                                }
                            }
                        }
                    }
                }
            }

            return 0;
        }

        PetscInt LinearElasticity::Hex8IsoparametricCE(PetscScalar dalpha, PetscScalar dbeta,PetscScalar* ce) {
            PetscErrorCode ierr=0;

            memset(ce, 0, sizeof(ce[0]) * 24 * 24);
            for (PetscInt i = 0; i < 24; i++) {
                for (PetscInt j = 0; j < 24; j++) {
                    ce[j + 24 * i] = dalpha*ME[j + 24 * i] + dbeta*KE[j + 24 * i] ;
                }
            }

            return ierr;
        }

        PetscScalar LinearElasticity::Dot(PetscScalar* v1, PetscScalar* v2, PetscInt l) {
            // Function that returns the dot product of v1 and v2,
            // which must have the same length l
            PetscScalar result = 0.0;
            for (PetscInt i = 0; i < l; i++) {
                result = result + v1[i] * v2[i];
            }
            return result;
        }

        void LinearElasticity::ShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta,
            PetscScalar* S) {
                // differentiatedShapeFunctions - Computes differentiated shape functions
                // At the point given by (xi, eta, zeta).
                // With respect to xi:
                S[0]=  0.125 * (1.0 - xi)*(1.0 - eta)*(1.0 - zeta);
                S[1]=  0.125 * (1.0 + xi)*(1.0 - eta)*(1.0 - zeta);
                S[2]=  0.125 * (1.0 + xi)*(1.0 + eta)*(1.0 - zeta);
                S[3]=  0.125 * (1.0 - xi)*(1.0 + eta)*(1.0 - zeta);
                S[4]=  0.125 * (1.0 - xi)*(1.0 - eta)*(1.0 + zeta);
                S[5]=  0.125 * (1.0 + xi)*(1.0 - eta)*(1.0 + zeta);
                S[6]=  0.125 * (1.0 + xi)*(1.0 + eta)*(1.0 + zeta);
                S[7]=  0.125 * (1.0 - xi)*(1.0 + eta)*(1.0 + zeta);

            }

            void LinearElasticity::ShapeFunctionsVector(PetscScalar* S,PetscScalar* vS) {

                vS[0]=  S[0];
                vS[1]=  0.0;
                vS[2]=  0.0;
                vS[3]=  S[1];
                vS[4]=  0.0;
                vS[5]=  0.0;
                vS[6]=  S[2];
                vS[7]=  0.0;
                vS[8]=  0.0;
                vS[9]=  S[3];
                vS[10]=  0.0;
                vS[11]=  0.0;
                vS[12]=  S[4];
                vS[13]=  0.0;
                vS[14]=  0.0;
                vS[15]=  S[5];
                vS[16]=  0.0;
                vS[17]=  0.0;
                vS[18]=  S[6];
                vS[19]=  0.0;
                vS[20]=  0.0;
                vS[21]=  S[7];
                vS[22]=  0.0;
                vS[23]=  0.0;

                vS[24]=  0.0;
                vS[25]=  S[0];
                vS[26]=  0.0;
                vS[27]=  0.0;
                vS[28]=  S[1];
                vS[29]=  0.0;
                vS[30]=  0.0;
                vS[31]=  S[2];
                vS[32]=  0.0;
                vS[33]=  0.0;
                vS[34]=  S[3];
                vS[35]=  0.0;
                vS[36]=  0.0;
                vS[37]=  S[4];
                vS[38]=  0.0;
                vS[39]=  0.0;
                vS[40]=  S[5];
                vS[41]=  0.0;
                vS[42]=  0.0;
                vS[43]=  S[6];
                vS[44]=  0.0;
                vS[45]=  0.0;
                vS[46]=  S[7];
                vS[47]=  0.0;

                vS[48]=  0.0;
                vS[49]=  0.0;
                vS[50]=  S[0];
                vS[51]=  0.0;
                vS[52]=  0.0;
                vS[53]=  S[1];
                vS[54]=  0.0;
                vS[55]=  0.0;
                vS[56]=  S[2];
                vS[57]=  0.0;
                vS[58]=  0.0;
                vS[59]=  S[3];
                vS[60]=  0.0;
                vS[61]=  0.0;
                vS[62]=  S[4];
                vS[63]=  0.0;
                vS[64]=  0.0;
                vS[65]=  S[5];
                vS[66]=  0.0;
                vS[67]=  0.0;
                vS[68]=  S[6];
                vS[69]=  0.0;
                vS[70]=  0.0;
                vS[71]=  S[7];


            }


            void LinearElasticity::DifferentiatedShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta,
                PetscScalar* dNdxi, PetscScalar* dNdeta, PetscScalar* dNdzeta) {
                    // differentiatedShapeFunctions - Computes differentiated shape functions
                    // At the point given by (xi, eta, zeta).
                    // With respect to xi:
                    dNdxi[0] = -0.125 * (1.0 - eta) * (1.0 - zeta);
                    dNdxi[1] = 0.125 * (1.0 - eta) * (1.0 - zeta);
                    dNdxi[2] = 0.125 * (1.0 + eta) * (1.0 - zeta);
                    dNdxi[3] = -0.125 * (1.0 + eta) * (1.0 - zeta);
                    dNdxi[4] = -0.125 * (1.0 - eta) * (1.0 + zeta);
                    dNdxi[5] = 0.125 * (1.0 - eta) * (1.0 + zeta);
                    dNdxi[6] = 0.125 * (1.0 + eta) * (1.0 + zeta);
                    dNdxi[7] = -0.125 * (1.0 + eta) * (1.0 + zeta);
                    // With respect to eta:
                    dNdeta[0] = -0.125 * (1.0 - xi) * (1.0 - zeta);
                    dNdeta[1] = -0.125 * (1.0 + xi) * (1.0 - zeta);
                    dNdeta[2] = 0.125 * (1.0 + xi) * (1.0 - zeta);
                    dNdeta[3] = 0.125 * (1.0 - xi) * (1.0 - zeta);
                    dNdeta[4] = -0.125 * (1.0 - xi) * (1.0 + zeta);
                    dNdeta[5] = -0.125 * (1.0 + xi) * (1.0 + zeta);
                    dNdeta[6] = 0.125 * (1.0 + xi) * (1.0 + zeta);
                    dNdeta[7] = 0.125 * (1.0 - xi) * (1.0 + zeta);
                    // With respect to zeta:
                    dNdzeta[0] = -0.125 * (1.0 - xi) * (1.0 - eta);
                    dNdzeta[1] = -0.125 * (1.0 + xi) * (1.0 - eta);
                    dNdzeta[2] = -0.125 * (1.0 + xi) * (1.0 + eta);
                    dNdzeta[3] = -0.125 * (1.0 - xi) * (1.0 + eta);
                    dNdzeta[4] = 0.125 * (1.0 - xi) * (1.0 - eta);
                    dNdzeta[5] = 0.125 * (1.0 + xi) * (1.0 - eta);
                    dNdzeta[6] = 0.125 * (1.0 + xi) * (1.0 + eta);
                    dNdzeta[7] = 0.125 * (1.0 - xi) * (1.0 + eta);
                }

                PetscScalar LinearElasticity::Inverse3M(PetscScalar J[][3], PetscScalar invJ[][3]) {
                    // inverse3M - Computes the inverse of a 3x3 matrix
                    PetscScalar detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2]) -
                    J[0][1] * (J[1][0] * J[2][2] - J[2][0] * J[1][2]) +
                    J[0][2] * (J[1][0] * J[2][1] - J[2][0] * J[1][1]);
                    invJ[0][0] = (J[1][1] * J[2][2] - J[2][1] * J[1][2]) / detJ;
                    invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) / detJ;
                    invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ;
                    invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) / detJ;
                    invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ;
                    invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) / detJ;
                    invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ;
                    invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) / detJ;
                    invJ[2][2] = (J[0][0] * J[1][1] - J[1][0] * J[0][1]) / detJ;
                    return detJ;
                }
