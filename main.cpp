
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "DarcyTest.h"

#include "TPZDarcyMaterial.h"

#include "TPZDarcyMaterial.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"


using namespace std;

const int SpaceHDiv          = 1; //Velocity in H(div) space
//const int SpaceContinuous  = 2; //Velocity [H1]ˆ2 space
//const int SpaceDiscontinuous = 3; //Velocity in H(Ph) - Ph: partition space

//const REAL Pi=M_PI;
const REAL visco        =  1.; // Coefficients: viscosity, permeability, symmetry factor
const REAL permeability =  1.; // Coefficients: viscosity, permeability, symmetry factor
const REAL theta        = -1.; // Coefficients: viscosity, permeability, symmetry factor

//const int matPoint =-5;//Material id of point
int dirichlet       = 0; // Conditions of problem boundary -> default Dirichlet left and right
int neumann         = 1; // Conditions of problem boundary -> default Dirichlet left and right
int penetration     = 2; // Conditions of problem boundary -> default Dirichlet left and right
int pointtype       = 5; // Conditions of problem boundary -> default Dirichlet left and right
int dirichletPress  = 6; // Conditions of problem boundary -> default Dirichlet left and right


// brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy, double r);


// brief Function to refine the geometric mesh
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);


// brief Function to create  the computational mesh for flux
TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder);


// brief Function to create  the computational mesh for pressure
TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);


// brief Function to create  the multi-physical computational mesh
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);


// Function to create interface between elements
TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);


// Function to add multiphysics interfaces
void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);



// definition of sol analytic
void sol_exact1(const TPZVec<REAL> & x, TPZVec<STATE>& f, TPZFMatrix<STATE> &dsol);

// Main function of the program:
int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    int h_level = 2;
    double hx   = 0.1; // Dimensions in X and Y of the domain
    double hy   = 1.0; // Dimensions in X and Y of the domain
    int nelx    = h_level; // Number of elements in X and Y
    int nely    = h_level; // Number of elements in X and Y
    int nx      = nelx+1; // Number of nodes in X and Y
    int ny      = nely+1; // Number of nodes in X and Y
    int pOrder  = 1; // Polynomial order of approximation
    
    
    DarcyTest * Test1 = new DarcyTest();
    Test1->Run(SpaceHDiv, pOrder, nx, ny, hx, hy,visco,permeability,theta);
    
    
  
//    std::cout << "FINISHED!" << std::endl;
//    
    return 0;
}






