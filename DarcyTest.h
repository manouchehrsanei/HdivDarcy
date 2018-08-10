//
//  DarcyTest.h
//  PZ
//
//  Created by Manouchehr on August 10, 2018.
//
//


#ifndef __PZ__DarcyTest__
#define __PZ__DarcyTest__

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
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"


using namespace std;
using namespace pzshape;

class DarcyTest{
private:
    
    int fdim; // Dimension of the problem
    int fmatID; // Material id
    
    // Material of contour
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;
    
    // Material of interface
    int fmatInterface;
    
    // Material of contour of interface
    int fmatIntBCbott;
    int fmatIntBCtop;
    int fmatIntBCleft;
    int fmatIntBCright;
    
    // Material of point
    int fmatPoint;
    
    // Boundary condition
    int fdirichlet;
    int fneumann;
    int fpenetration;
    int fpointtype;
    int fdirichletvar;
    
    
    int fquadmat1; //Bottom of the square
    int fquadmat2; //Top of the square
    int fquadmat3; //Material of interface
    
    STATE fviscosity;
    STATE fpermeability;
    STATE ftheta;
    
    
public:
    
    bool fisH1;
    
    
    DarcyTest();
    
    ~DarcyTest();
    
    void Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE permeability, STATE theta);
    
    /*  Geometry mesh */
    TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy);
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    
    /* Computational mesh */
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
    
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE permeability, STATE theta);
    
    
    // Exact solution
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    // Right side of the equation
    static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    // static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    
    
};


#endif
