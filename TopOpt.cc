#include <TopOpt.h>
#include <cmath>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
 Updated: June 2019, Niels Aage
 Copyright (C) 2013-2019,

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

TopOpt::TopOpt(PetscInt nconstraints) {

    m = nconstraints;
    Init();
}

TopOpt::TopOpt() {

    m = 1;
    Init();
}

void TopOpt::Init() { // Dummy constructor

    x        = NULL;
    xPhys    = NULL;
    dfdx     = NULL;
    dgdx     = NULL;
    gx       = NULL;
    da_nodes = NULL;
    da_elem  = NULL;

    xo1 = NULL;
    xo2 = NULL;
    U   = NULL;
    L   = NULL;

    SetUp();
}

TopOpt::~TopOpt() {

    // Delete vectors
    if (x != NULL) {VecDestroy(&x);}
    if (xTilde != NULL) {VecDestroy(&xTilde);}
    if (xPhys != NULL) {VecDestroy(&xPhys);}
    if (dfdx != NULL) {VecDestroy(&dfdx);}
    if (dgdx != NULL) {VecDestroyVecs(m, &dgdx);}
    if (xold != NULL) {VecDestroy(&xold);}
    if (xmin != NULL) {VecDestroy(&xmin);}
    if (xmax != NULL) {VecDestroy(&xmax);}

    if (da_nodes != NULL) {DMDestroy(&(da_nodes));}
    if (da_elem != NULL) {DMDestroy(&(da_elem));}

    // Delete constraints
    if (gx != NULL) {delete[] gx;}

    // mma restart method
    if (xo1 != NULL) {VecDestroy(&xo1);}
    if (xo2 != NULL) {VecDestroy(&xo2);}
    if (L != NULL) {VecDestroy(&L);}
    if (U != NULL) {VecDestroy(&U);}
}

// PetscErrorCode TopOpt::SetUp(Vec CRAPPY_VEC){
PetscErrorCode TopOpt::SetUp() {
    PetscErrorCode ierr;

      PetscBool flg;
	boundaries = 0;
    PetscOptionsGetInt(NULL, NULL, "-boundaries", &boundaries, &flg);
    
    switch (boundaries){
        case 0: // MBB type of boundary conditions 

    nxyz[0] = 97; // 129;
    nxyz[1] = 33; // 65;
    nxyz[2] = 33; // 65;

	xc[0]   = 0.0;
    xc[1]   = 60.0;
    xc[2]   = 0.0;
    xc[3]   = 20.0;
    xc[4]   = 0.0;
    xc[5]   = 20.0;
	 loadloc = 0;
        break;
        case 1: // Cantilever type of boundary conditions


    nxyz[0] = 65 ; // 129;
    nxyz[1] = 33; // 65;
    nxyz[2] = 33; // 65;
	 
	 xc[0]   = 0.0;
    xc[1]   = 40.0;
    xc[2]   = 0.0;
    xc[3]   = 20.0;
    xc[4]   = 0.0;
    xc[5]   = 20.0;
	 loadloc = 1;
        break;
    }
    nu      = 0.3;
    nlvls   = 4;

    // SET DEFAULTS for optimization problems
    volfrac = 0.3; //volume fraction constraint
    V0 = volfrac; // initial volume fraction
    maxItr  = 10;
    rmin    = 1.5;
    penalK   = 3.0; // Use 3 for simp, 5 for RAMP
    penalM   = 1.0;
    Emax    = 1E3;
    Rmax    = 4.16666666667e-4; //4.166666666667e-4; // Total structure has a weight of 10kg
    filter  = 1; // 0=sens,1=dens,2=PDE - other val == no filtering
	 sInt = 1; //0: linear, 1: SIMP, 2: RAMP
	 mInt = 0;			// 0: linear, 1: SIMP, 2: RAMP
    Xmin    = 0.0;
    Xmax    = 1.0;
    movlim  = 0.2;
    //restart = PETSC_FALSE;
    Rf = 1e-3;
    Ef = 1e-9;
    Emin = Ef*Emax;
    Rmin = Rf*Rmax;

    ti  = 1; // Time integration (0: backward Euler, 1: Newmark )
    t0  = 0.0; // Initial time
    t1  = 50 ; // End time

    ncp = 2; // NUmber of checkpoints
    tsc = 25; // number of timesteps per checkpoint

    dalpha = 0.0001; // Damping alpha
    dbeta = 0.0001; // damping beta


	//Objective function
	objective = 0;
	minmax = -1; // default minimize 
	boundaries = 0; // 0: MBB, 1: cantilever 
	optimizationRegion = 0; // Default, all elements 

    // Projection filter
    projectionFilter = PETSC_TRUE;
    beta             = 0.1;
    betaFinal        = 48;
    eta              = 0.15;

	 // SOAR settings
	 reduction  = PETSC_FALSE;
	 basis = 5;
	 omega0 = 1.0; // HAS UNITS Hz!!! 

    ierr = SetUpMESH();     CHKERRQ(ierr);
    ierr = SetUpOPT();      CHKERRQ(ierr);

    return (ierr);
}

PetscErrorCode TopOpt::SetUpMESH() {

    PetscErrorCode ierr;

    // Read input from arguments
    PetscBool flg;

    // Physics parameters

    // Write parameters for the physics _ OWNED BY TOPOPT
    // Optimization paramteres
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# M E S H   S E T T I N G S                              \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");
    PetscOptionsGetInt(NULL, NULL, "-nx", &(nxyz[0]), &flg);
    PetscOptionsGetInt(NULL, NULL, "-ny", &(nxyz[1]), &flg);
    PetscOptionsGetInt(NULL, NULL, "-nz", &(nxyz[2]), &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Number of nodes: (-nx,-ny,-nz):        (%i,%i,%i) \n", nxyz[0], nxyz[1], nxyz[2]);
    PetscPrintf(PETSC_COMM_WORLD, "# Number of degree of freedom:           %i \n", 3 * nxyz[0] * nxyz[1] * nxyz[2]);
    PetscPrintf(PETSC_COMM_WORLD, "# Number of elements:                    (%i,%i,%i) \n", nxyz[0] - 1, nxyz[1] - 1,nxyz[2] - 1);
    PetscOptionsGetReal(NULL, NULL, "-xcmin", &(xc[0]), &flg);
    PetscOptionsGetReal(NULL, NULL, "-xcmax", &(xc[1]), &flg);
    PetscOptionsGetReal(NULL, NULL, "-ycmin", &(xc[2]), &flg);
    PetscOptionsGetReal(NULL, NULL, "-ycmax", &(xc[3]), &flg);
    PetscOptionsGetReal(NULL, NULL, "-zcmin", &(xc[4]), &flg);
    PetscOptionsGetReal(NULL, NULL, "-zcmax", &(xc[5]), &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Dimensions: (-xcmin,-xcmax,..,-zcmax): (%.2f,%.2f,%.2f)\n", xc[1] - xc[0], xc[3] - xc[2],xc[5] - xc[4]);
    PetscPrintf(PETSC_COMM_WORLD, "# elemnt size: (dx, dy, dz)            : (%.2e,%.2e,%.2e)\n", (xc[1] - xc[0])/(nxyz[0] - 1), (xc[3] - xc[2])/(nxyz[1] - 1),(xc[5] - xc[4])/(nxyz[2] - 1));
    PetscOptionsGetInt(NULL, NULL, "-nlvls", &nlvls,&flg); // NEEDS THIS TO CHECK IF MESH IS OK BEFORE PROCEEDING !!!!
    PetscPrintf(PETSC_COMM_WORLD, "# -nlvls: %i\n", nlvls);
    PetscPrintf(PETSC_COMM_WORLD, "#\n");


    // Check if the mesh supports the chosen number of MG levels
    PetscScalar divisor = PetscPowScalar(2.0, (PetscScalar)nlvls - 1.0);
    // x - dir
    if (std::floor((PetscScalar)(nxyz[0] - 1) / divisor) != (nxyz[0] - 1.0) / ((PetscInt)divisor)) {
        PetscPrintf(PETSC_COMM_WORLD, "MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
        PetscPrintf(PETSC_COMM_WORLD, "X - number of nodes %i is cannot be halfened %i times\n", nxyz[0], nlvls - 1);
        exit(0);
    }
    // y - dir
    if (std::floor((PetscScalar)(nxyz[1] - 1) / divisor) != (nxyz[1] - 1.0) / ((PetscInt)divisor)) {
        PetscPrintf(PETSC_COMM_WORLD, "MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
        PetscPrintf(PETSC_COMM_WORLD, "Y - number of nodes %i is cannot be halfened %i times\n", nxyz[1], nlvls - 1);
        exit(0);
    }
    // z - dir
    if (std::floor((PetscScalar)(nxyz[2] - 1) / divisor) != (nxyz[2] - 1.0) / ((PetscInt)divisor)) {
        PetscPrintf(PETSC_COMM_WORLD, "MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
        PetscPrintf(PETSC_COMM_WORLD, "Z - number of nodes %i is cannot be halfened %i times\n", nxyz[2], nlvls - 1);
        exit(0);
    }


    // Start setting up the FE problem
    // Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED,
    // DMDA_BOUNDARY_PERIODIC
    DMBoundaryType bx = DM_BOUNDARY_NONE;
    DMBoundaryType by = DM_BOUNDARY_NONE;
    DMBoundaryType bz = DM_BOUNDARY_NONE;

    // Stencil type - box since this is closest to FEM (i.e. STAR is FV/FD)
    DMDAStencilType stype = DMDA_STENCIL_BOX;

    // Discretization: nodes:
    // For standard FE - number must be odd
    // FOr periodic: Number must be even
    PetscInt nx = nxyz[0];
    PetscInt ny = nxyz[1];
    PetscInt nz = nxyz[2];

    // number of nodal dofs: Nodal design variable - NOT REALLY NEEDED
    PetscInt numnodaldof = 1;

    // Stencil width: each node connects to a box around it - linear elements
    PetscInt stencilwidth = 1;

    // Coordinates and element sizes: note that dx,dy,dz are half the element size
    PetscReal xmin = xc[0], xmax = xc[1], ymin = xc[2], ymax = xc[3], zmin = xc[4], zmax = xc[5];
    dx = (xc[1] - xc[0]) / (PetscScalar(nxyz[0] - 1));
    dy = (xc[3] - xc[2]) / (PetscScalar(nxyz[1] - 1));
    dz = (xc[5] - xc[4]) / (PetscScalar(nxyz[2] - 1));

    // Create the nodal mesh
    ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, stype, nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        numnodaldof, stencilwidth, 0, 0, 0, &(da_nodes));
    CHKERRQ(ierr);

    // Initialize
    DMSetFromOptions(da_nodes);
    DMSetUp(da_nodes);

    // Set the coordinates
    ierr = DMDASetUniformCoordinates(da_nodes, xmin, xmax, ymin, ymax, zmin, zmax);
    CHKERRQ(ierr);

    // Set the element type to Q1: Otherwise calls to GetElements will change to
    // P1 ! STILL DOESN*T WORK !!!!
    ierr = DMDASetElementType(da_nodes, DMDA_ELEMENT_Q1);
    CHKERRQ(ierr);

    // Create the element mesh: NOTE THIS DOES NOT INCLUDE THE FILTER !!!
    // find the geometric partitioning of the nodal mesh, so the element mesh will
    // coincide with the nodal mesh
    PetscInt md, nd, pd;
    ierr = DMDAGetInfo(da_nodes, NULL, NULL, NULL, NULL, &md, &nd, &pd, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);

    // vectors with partitioning information
    PetscInt* Lx = new PetscInt[md];
    PetscInt* Ly = new PetscInt[nd];
    PetscInt* Lz = new PetscInt[pd];

    // get number of nodes for each partition
    const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
    ierr = DMDAGetOwnershipRanges(da_nodes, &LxCorrect, &LyCorrect, &LzCorrect);
    CHKERRQ(ierr);

    // subtract one from the lower left corner.
    for (int i = 0; i < md; i++) {
        Lx[i] = LxCorrect[i];
        if (i == 0) {
            Lx[i] = Lx[i] - 1;
        }
    }
    for (int i = 0; i < nd; i++) {
        Ly[i] = LyCorrect[i];
        if (i == 0) {
            Ly[i] = Ly[i] - 1;
        }
    }
    for (int i = 0; i < pd; i++) {
        Lz[i] = LzCorrect[i];
        if (i == 0) {
            Lz[i] = Lz[i] - 1;
        }
    }

    // Create the element grid: NOTE CONNECTIVITY IS 0
    PetscInt conn = 0;
    ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, stype, nx - 1, ny - 1, nz - 1, md, nd, pd, 1, conn, Lx, Ly, Lz,
                        &(da_elem));
    CHKERRQ(ierr);

    // Initialize
    DMSetFromOptions(da_elem);
    DMSetUp(da_elem);

    // Set element center coordinates
    ierr = DMDASetUniformCoordinates(da_elem, xmin + dx / 2.0, xmax - dx / 2.0, ymin + dy / 2.0, ymax - dy / 2.0,
                                     zmin + dz / 2.0, zmax - dz / 2.0);
    CHKERRQ(ierr);

    // Clean up
    delete[] Lx;
    delete[] Ly;
    delete[] Lz;

    return (ierr);
}

PetscErrorCode TopOpt::SetUpOPT() {

    PetscErrorCode ierr;

    // ierr = VecDuplicate(CRAPPY_VEC,&xPhys); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da_elem, &xPhys);
    CHKERRQ(ierr);
    // Total number of design variables
    VecGetSize(xPhys, &n);

    PetscBool flg;

    // Optimization paramteres
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# B O U N D A R Y  C O N D I T I O N S                   \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");
    PetscOptionsGetInt(NULL, NULL, "-boundaries", &boundaries, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Boundary conditions (-boundaries): \n");
	PetscPrintf(PETSC_COMM_WORLD,"#   0: MBB TYPE        :  "); if (boundaries == 0){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   1: CANTILEVER TYPE :  "); if (boundaries == 1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
   PetscPrintf(PETSC_COMM_WORLD, "# \n");
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# L O A D  L O C A T I O N                               \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");
    PetscOptionsGetInt(NULL, NULL, "-loadloc", &loadloc, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Load location preset (-loadloc): \n");
	PetscPrintf(PETSC_COMM_WORLD,"#   0: MBB TYPE        :  "); if (loadloc == 0){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   1: CANTILEVER TYPE :  "); if (loadloc == 1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# O P T I M I Z A T I O N   S E T T I N G S              \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Problem size: n= %i, m= %i\n", n, m);
    PetscOptionsGetInt(NULL, NULL, "-filter", &filter, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Filter (-filter) (0=sens., 1=dens, 2=PDE)   : %i  \n", filter);
    PetscOptionsGetReal(NULL, NULL, "-rmin", &rmin, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Filter radius (-rmin)                       : %f \n", rmin);
    PetscOptionsGetBool(NULL, NULL, "-projectionFilter", &projectionFilter, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Projection filter (-projectionFilter) (0/1) : %i \n", projectionFilter);
    PetscOptionsGetReal(NULL, NULL, "-beta", &beta, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Initial beta (-beta)                        : %f \n", beta);
    PetscOptionsGetReal(NULL, NULL, "-betaFinal", &betaFinal, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Final Beta (-betaFinal)                     : %f \n", betaFinal);
    PetscOptionsGetReal(NULL, NULL, "-eta", &eta, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Projection threshold (-eta)                 : %f \n", eta);
    PetscOptionsGetReal(NULL, NULL, "-volfrac", &volfrac, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Volume fraction (-volfrac)                  : %f \n", volfrac);
    PetscOptionsGetReal(NULL, NULL, "-Xmax", &Xmax, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Xmax (-Xmax)                                : %f \n", Xmax);
    PetscOptionsGetReal(NULL, NULL, "-Xmin", &Xmin, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Xmin (-Xmin)                                : %f \n", Xmin);
    PetscOptionsGetInt(NULL, NULL, "-maxItr", &maxItr, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Iter max (-maxItr)                          : %i \n", maxItr);
    PetscOptionsGetReal(NULL, NULL, "-movlim", &movlim, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Move limit (-movlim)                        : %f \n", movlim);
    PetscOptionsGetInt(NULL, NULL, "-optimizationRegion", &optimizationRegion, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Optimization Region preset (-optimizationRegion): \n");
	PetscPrintf(PETSC_COMM_WORLD,"#   0: Vicinity of -loadloc 0  :  "); if (optimizationRegion == 0){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   1: Vicinity of -loadloc 1  :  "); if (optimizationRegion == 1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   2: All elements            :  "); if (optimizationRegion == 2){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   3: 5x5x5 box inside domain :  "); if (optimizationRegion == 3){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
    PetscPrintf(PETSC_COMM_WORLD, "# \n");
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# M A T E R I A L S  \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");

    PetscOptionsGetReal(NULL, NULL, "-V0", &V0, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Initial volume fraction (-V0)               : %f \n", V0);
    PetscOptionsGetReal(NULL, NULL, "-nu", &nu, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Poisson's ratio       (-nu)                 : %f \n", nu);
    PetscOptionsGetReal(NULL, NULL, "-penalM", &penalM, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Mass penalization (-penalM)                 : %f \n", penalM);
    PetscOptionsGetReal(NULL, NULL, "-penalK", &penalK, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Stiffness penalization (-penalK)            : %f \n", penalK);
    PetscOptionsGetReal(NULL, NULL, "-Rf", &Rf, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Contrast, mass (-Rf)                        : %e \n", Rf);
    PetscOptionsGetReal(NULL, NULL, "-Ef", &Ef, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Contrast, Stiffness (-Ef)                   : %e \n", Ef);
    PetscOptionsGetReal(NULL, NULL, "-Emax", &Emax, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Solid, stiffness (-Emax)                    : %e \n", Emax);
    PetscOptionsGetReal(NULL, NULL, "-Rmax", &Rmax, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Solid, Mass (-Rmax)                         : %e \n", Rmax);
    Emin = Ef*Emax;
    PetscPrintf(PETSC_COMM_WORLD, "# Void, stiffness (Ef x Emax)                 : %e \n", Emin);
    Rmin = Rf*Rmax;
    PetscPrintf(PETSC_COMM_WORLD, "# Void, Mass (Rf x Rmax)                      : %e \n", Rmin);
    PetscOptionsGetReal(NULL, NULL, "-dalph", &dalpha, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Damping, alpha (-dalpha)                    : %e \n", dalpha);
    PetscOptionsGetReal(NULL, NULL, "-dbeta", &dbeta, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Damping, beta (-dbeta)                      : %e \n", dbeta);
    PetscOptionsGetInt(NULL, NULL, "-massInterpolation", &mInt, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# MassInterpolation scheme: (-mInt) (0=Linear, 1=SIMP, 2=RAMP)   : %i  \n", mInt);
    PetscOptionsGetInt(NULL, NULL, "-sInt", &sInt, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# StiffnessInterpolation scheme: (-sInt) (0=Linear., 1=SIMP, 2=RAMP)   : %i  \n", sInt);
    PetscPrintf(PETSC_COMM_WORLD, "# \n");
    PetscPrintf(PETSC_COMM_WORLD, "#########################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# T E M P O R A L P R O P E R T I E S                    \n");
    PetscPrintf(PETSC_COMM_WORLD, "# -------------------------------------------------------\n");
    PetscOptionsGetInt(NULL, NULL, "-ti", &ti, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Timeintegration scheme (-ti)  (0=BE, 1=NM)  : %i  \n", ti);
    PetscOptionsGetReal(NULL, NULL, "-t0", &t0, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Initial time                                : %f  \n", t0);
    PetscOptionsGetReal(NULL, NULL, "-t1", &t1, &flg);
    PetscPrintf(PETSC_COMM_WORLD, "# Final time             (-t1)                : %f  \n", t1);
    PetscOptionsGetInt(NULL, NULL, "-ncp", &ncp, &flg);
    PetscOptionsGetInt(NULL, NULL, "-tsc", &tsc, &flg);
    Ts = ncp*tsc+1;                        // Number of timesteps
    PetscPrintf(PETSC_COMM_WORLD, "# Number of checkpoints  (-ncp)               : %i  \n", ncp);
    PetscPrintf(PETSC_COMM_WORLD, "# Number of timesteps/checkpoint (-tsc)       : %i  \n", tsc);
    PetscPrintf(PETSC_COMM_WORLD, "# Number of timesteps                         : %i  \n", Ts);
    PetscScalar dt = t1/(Ts-1);
    PetscPrintf(PETSC_COMM_WORLD, "# Timestep size (t1 / (Ts-1))                 : %f  \n", dt);
    PetscPrintf(PETSC_COMM_WORLD, "# \n");
   PetscPrintf(PETSC_COMM_WORLD,"#########################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# O B J E C T I V E  F U N C T I O N\n");
	PetscPrintf(PETSC_COMM_WORLD,"#--------------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD,"#--------------------------------------------------------\n");
    PetscOptionsGetInt(NULL, NULL, "-objective", &objective, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Objective function (-objective): \n");
	PetscPrintf(PETSC_COMM_WORLD,"#   0: Squared displacements :  "); if (objective == 0){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   1: Squared velocities    :  "); if (objective == 1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   2: Squared accelerations :  "); if (objective == 2){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   3: Potential energy      :  "); if (objective == 3){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   4: Kinetic energy        :  "); if (objective == 4){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   5: Total energy          :  "); if (objective == 5){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   6: Dissipated Energy     :  "); if (objective == 6){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
   PetscPrintf(PETSC_COMM_WORLD, "# \n");
    PetscOptionsGetInt(NULL, NULL, "-minmax", &minmax, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Mininimize/maximise (-minmax) : \n");
	PetscPrintf(PETSC_COMM_WORLD,"#   Minimize:  "); ; if (minmax == -1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
	PetscPrintf(PETSC_COMM_WORLD,"#   Maximize:  "); ; if (minmax == 1){PetscPrintf(PETSC_COMM_WORLD,"1\n");} else {PetscPrintf(PETSC_COMM_WORLD,"0\n");}
   PetscPrintf(PETSC_COMM_WORLD,"#########################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# R E D U C T I O N                                      \n");
	PetscPrintf(PETSC_COMM_WORLD,"#--------------------------------------------------------\n");
    PetscOptionsGetBool(NULL, NULL, "-reduction", &reduction, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Reduction (1/0) (-reduction): %i\n",reduction);
    PetscOptionsGetInt(NULL, NULL, "-basis", &basis, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Number of basis vectors (-basis): %i\n",basis);
    PetscOptionsGetScalar(NULL, NULL, "-omega0", &omega0, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# Omega0 (-omega0) *IMPORTANT*: %f [Hz]\n",omega0);
	PetscPrintf(PETSC_COMM_WORLD,"# Omega0 (-omega0) *IMPORTANT*: %f [rad]\n",2*pi*omega0);
	PetscPrintf(PETSC_COMM_WORLD,"# Frequencyoffset for reduction 0.1 x omega0 : %f [Hz] \n",0.1*omega0);
   PetscPrintf(PETSC_COMM_WORLD, "# \n");
   PetscPrintf(PETSC_COMM_WORLD,"#########################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# O U T P U T S                                          \n");
	PetscPrintf(PETSC_COMM_WORLD,"#--------------------------------------------------------\n");
	outputPhysics0 = PETSC_TRUE; // default
	outputPhysics1 = PETSC_TRUE; // default
	outputDesign = PETSC_TRUE; // default
PetscOptionsGetBool(NULL, NULL, "-outputPhysics0", &outputPhysics0, &flg);
PetscOptionsGetBool(NULL, NULL, "-outputPhysics1", &outputPhysics1, &flg);
PetscOptionsGetBool(NULL, NULL, "-outputDesign", &outputDesign, &flg);
	PetscPrintf(PETSC_COMM_WORLD,"# (-outputPhysics0): %d \n",outputPhysics0); //state of initial structure incl design 
	PetscPrintf(PETSC_COMM_WORLD,"# (-outputPhysics1): %d \n",outputPhysics1); //State of optimized structure incl design
	PetscPrintf(PETSC_COMM_WORLD,"# (-outputDesign)  : %d \n",outputDesign); //design iterations 
   PetscPrintf(PETSC_COMM_WORLD, "# \n");
   PetscPrintf(PETSC_COMM_WORLD,"#########################################################\n");
   PetscPrintf(PETSC_COMM_WORLD, "# \n");
  
  
  // Allocate after input
    gx = new PetscScalar[m];
    if (filter == 0) {
        Xmin = 0.001; // Prevent division by zero in filter
    }

    // Allocate the optimization vectors
    ierr = VecDuplicate(xPhys, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(xPhys, &xTilde);
    CHKERRQ(ierr);

    VecSet(x, V0);      // Initialize to volfrac !
    VecSet(xTilde, V0); // Initialize to volfrac !
    VecSet(xPhys, V0);  // Initialize to volfrac


    // Sensitivity vectors
    ierr = VecDuplicate(x, &dfdx);
    CHKERRQ(ierr);
    ierr = VecDuplicateVecs(x, m, &dgdx);
    CHKERRQ(ierr);

    // Bounds and
    VecDuplicate(x, &xmin);
    VecDuplicate(x, &xmax);
    VecDuplicate(x, &xold);
    VecSet(xold, V0);

    return (ierr);
}

PetscErrorCode TopOpt::AllocateMMAwithRestart(PetscInt* itr, MMA** mma) {
    PetscErrorCode ierr = 0;

    // Set MMA parameters (for multiple load cases)
    PetscScalar aMMA[m];
    PetscScalar cMMA[m];
    PetscScalar dMMA[m];
    for (PetscInt i = 0; i < m; i++) {
        aMMA[i] = 0.0;
        dMMA[i] = 0.0;
        cMMA[i] = 1000.0;
    }

    // Check if restart is desired
    restart                  = PETSC_TRUE;  // DEFAULT USES RESTART
    flip                     = PETSC_TRUE;  // BOOL to ensure that two dump streams are kept
    PetscBool onlyLoadDesign = PETSC_FALSE; // Default restarts everythingi
    PetscInt loadDesign = 2; // 0: load x, 1: load xT, 2: load xP

    // Get inputs
    PetscBool flg;
    char      filenameChar[PETSC_MAX_PATH_LEN];
    PetscOptionsGetBool(NULL, NULL, "-restart", &restart, &flg);
    PetscOptionsGetBool(NULL, NULL, "-onlyLoadDesign", &onlyLoadDesign, &flg);
    PetscOptionsGetInt(NULL, NULL, "-loadDesign", &loadDesign, &flg);

    if (restart) {
        ierr = VecDuplicate(x, &xo1);
        CHKERRQ(ierr);
        ierr = VecDuplicate(x, &xo2);
        CHKERRQ(ierr);
        ierr = VecDuplicate(x, &U);
        CHKERRQ(ierr);
        ierr = VecDuplicate(x, &L);
        CHKERRQ(ierr);
    }

    // Determine the right place to write the new restart files
    std::string filenameWorkdir = "./";
    PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {

        filenameWorkdir = "";
        filenameWorkdir.append(filenameChar);
    }
	
    filename00    = filenameWorkdir;
    filename00Itr = filenameWorkdir;
    filename01    = filenameWorkdir;
    filename01Itr = filenameWorkdir;

    filename00.append("/Restart00.dat");
    filename00Itr.append("/Restart00_itr_f0.dat");
    filename01.append("/Restart01.dat");
    filename01Itr.append("/Restart01_itr_f0.dat");

    // Where to read the restart point from
    std::string restartFileVec = ""; // NO RESTART FILE !!!!!
    std::string restartFileItr = ""; // NO RESTART FILE !!!!!

    PetscOptionsGetString(NULL, NULL, "-restartFileVec", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {
        restartFileVec.append(filenameChar);
    }
    PetscOptionsGetString(NULL, NULL, "-restartFileItr", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {
        restartFileItr.append(filenameChar);
    }

    // Which solution to use for restarting
    PetscPrintf(PETSC_COMM_WORLD, "##############################################################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Continue from previous iteration (-restart): %i \n", restart);
    PetscPrintf(PETSC_COMM_WORLD, "# Restart file (-restartFileVec): %s \n", restartFileVec.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "# Restart file (-restartFileItr): %s \n", restartFileItr.c_str());
    PetscPrintf(PETSC_COMM_WORLD,
                "# New restart files are written to (-workdir): %s "
                "(Restart0x.dat and Restart0x_itr_f0.dat) \n",
                filenameWorkdir.c_str());

    // Check if files exist:
    PetscBool vecFile = fexists(restartFileVec);
    if (!vecFile) {
        PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileVec.c_str());
    }
    PetscBool itrFile = fexists(restartFileItr);
    if (!itrFile) {
        PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileItr.c_str());
    }

    // Read from restart point
    PetscScalar RHSnorm;
    PetscInt nGlobalDesignVar;
    VecGetSize(x,
               &nGlobalDesignVar); // ASSUMES THAT SIZE IS ALWAYS MATCHED TO CURRENT MESH
    if (restart && vecFile && itrFile) {


        PetscViewer view;
        // Open the data files
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, restartFileVec.c_str(), FILE_MODE_READ, &view);

        VecLoad(x, view);
		  VecSum(x,&RHSnorm); PetscPrintf(PETSC_COMM_WORLD,"# sum(x) LOADED: %8f\n",RHSnorm);
        VecLoad(xPhys, view);
		  VecSum(xPhys,&RHSnorm); PetscPrintf(PETSC_COMM_WORLD,"# sum(xPhys) LOADED: %16f\n",RHSnorm-796215.021281);
        VecLoad(xo1, view);
        VecLoad(xo2, view);
        VecLoad(U, view);
        VecLoad(L, view);
        PetscViewerDestroy(&view);

	 // Read iteration and fscale
        std::fstream itrfile(restartFileItr.c_str(), std::ios_base::in);
//        itrfile >> itr[0];
//        itrfile >> fscale;

        // Choose if restart is full or just an initial design guess
        if (onlyLoadDesign) {	

				switch (loadDesign){	
				case 0: 
					// dont do anything, x is loaded in lines below. 
					PetscPrintf(PETSC_COMM_WORLD,"# x was loaded from resart file loaded, as unsed as initial x.\n");
				break;
				case 1:
					// Is xP even available?
					PetscPrintf(PETSC_COMM_WORLD,"# xT was (not?) loaded from resart file loaded, and used as initial x.\n");
				break;
				case 2:	
					PetscPrintf(PETSC_COMM_WORLD,"# xPhys was loads from resart file loaded, and used as initial x.\n");
					VecCopy(xPhys,x);
				break;
				}
				
            PetscPrintf(PETSC_COMM_WORLD, "# Loading design from file: %s \n", restartFileVec.c_str());
            *mma = new MMA(nGlobalDesignVar, m, x, aMMA, cMMA, dMMA);

        } else {
		  VecCopy(x,xold);
        itrfile >> itr[0];
        itrfile >> fscale;
            PetscPrintf(PETSC_COMM_WORLD, "# Continue optimization from file: %s \n", restartFileVec.c_str());
            *mma = new MMA(nGlobalDesignVar, m, *itr, xo1, xo2, U, L, aMMA, cMMA, dMMA);
   
		PetscPrintf(PETSC_COMM_WORLD,"# written fscale: %1.16f\n",fscale);
		PetscPrintf(PETSC_COMM_WORLD,"# written itr[0]: %i\n",itr[0]);

        }

        PetscPrintf(PETSC_COMM_WORLD, "# Successful restart from file: %s and %s \n", restartFileVec.c_str(),
                    restartFileItr.c_str());
    } else {
        *mma = new MMA(nGlobalDesignVar, m, x, aMMA, cMMA, dMMA);
    }

	PetscPrintf(PETSC_COMM_WORLD,"# Finished allocating MMAwithrestart\n");
    return ierr;
}

PetscErrorCode TopOpt::WriteRestartFiles(PetscInt* itr, MMA* mma) {

    PetscErrorCode ierr = 0;
    // Only dump data if correct allocater has been used
    if (!restart) {
        return -1;
    }

    // Get restart vectors
    mma->Restart(xo1, xo2, U, L);
 
 // Choose previous set of restart files
    if (flip) {
        flip = PETSC_FALSE;
    } else {
        flip = PETSC_TRUE;
    }

    // Write file with iteration number of f0 scaling
    // and a file with the MMA-required vectors, in the following order:
    // : x,xPhys,xold1,xold2,U,L
    PetscViewer view;         // vectors
    PetscViewer restartItrF0; // scalars

    PetscViewerCreate(PETSC_COMM_WORLD, &restartItrF0);
    PetscViewerSetType(restartItrF0, PETSCVIEWERASCII);
    PetscViewerFileSetMode(restartItrF0, FILE_MODE_WRITE);

    // Open viewers for writing
    if (!flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename00.c_str(), FILE_MODE_WRITE, &view);
        PetscViewerFileSetName(restartItrF0, filename00Itr.c_str());
    } else if (flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename01.c_str(), FILE_MODE_WRITE, &view);
        PetscViewerFileSetName(restartItrF0, filename01Itr.c_str());
    }

    // Write iteration and fscale
    PetscViewerASCIIPrintf(restartItrF0, "%d ", itr[0]);
    PetscViewerASCIIPrintf(restartItrF0, "%1.16e", fscale);
    PetscViewerASCIIPrintf(restartItrF0, "\n");

    // Write vectors
    VecView(x, view); // the design variables
    VecView(xPhys, view);
    VecView(xo1, view);
    VecView(xo2, view);
    VecView(U, view);
    VecView(L, view);

    // Clean up
    PetscViewerDestroy(&view);
    PetscViewerDestroy(&restartItrF0);

    // PetscPrintf(PETSC_COMM_WORLD,"DONE WRITING DATA\n");
    return ierr;
}
