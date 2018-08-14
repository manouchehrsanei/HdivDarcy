
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

int m_dim           =  2; // Dimension of the problem
STATE fviscosity    =  1.;
STATE fpermeability =  1.;
STATE ftheta        = -1.;

// Geometry

int pOrder          =  2; // Polynomial order of approximation
int h_level         =  8;
double hx           =  1.0; // Dimensions in X and Y of the domain
double hy           =  1.0; // Dimensions in X and Y of the domain
int nelx            =  h_level; // Number of elements in X and Y
int nely            =  h_level; // Number of elements in X and Y
int nx              =  nelx + 1; // Number of nodes in X and Y
int ny              =  nely + 1; // Number of nodes in X and Y

int m_matID         =  1; // Material id

// Boundary condition
int dirichlet       =  0;
int neumann         =  1;
int penetration     =  2;
int pointtype       =  5;
int dirichletPress  =  6;

// Boundary condition
int m_dirichlet      = 0;
int m_neumann        = 1;
int m_penetration    = 2;
int m_pointtype      = 3;
int m_dirichletvar   = 4;


// Material of contour
int m_matBCbott      = -1;
int m_matBCtop       = -2;
int m_matBCleft      = -3;
int m_matBCright     = -4;

// Material of interface
int m_matInterface   = 4;

// Material of contour of interface
int m_matIntBCbott   = -11;
int m_matIntBCtop    = -12;
int m_matIntBCleft   = -13;
int m_matIntBCright  = -14;

// Material of point
int m_matPoint       = -5;

//The square
int m_quadmat1        =  1;
int m_quadmat2        =  2;
int m_quadmat3        =  3;


// brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh(int nelx, int nely, double hx, double hy);


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

// Exact solution
static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);

// Right side of the equation
static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);



// Main function of the program:
int main()
{
    
    TPZMaterial::gBigNumber = 1.e16;
    HDivPiola = 1;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    // Create geomesh
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy);
    
#ifdef PZDEBUG
//    std::ofstream meshfile("GeoMesh.txt");
//    gmesh->Print(meshfile);
//    
//    std::ofstream meshvtk("GeoMesh.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, meshvtk,true);
#endif
    
    // Create computational mesh
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); // Computational mesh for velocity
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); // Computational mesh for pressure
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); // Multiphysics mesh
    
#ifdef PZDEBUG
//    std::ofstream comVfile("CompMesh_V.txt");
//    cmesh_v->Print(comVfile);
//    
//    std::ofstream comPfile("CompMesh_P.txt");
//    cmesh_p->Print(comPfile);
//    
//    std::ofstream comMfile("CompMesh_M.txt");
//    cmesh_m->Print(comMfile);
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; // Prevents of renumbering of the equations (As the same of Oden's result)
    TPZAnalysis analysis(cmesh_m, optimizeBandwidth); // Create analysis
    
    TPZSkylineStructMatrix matskl(cmesh_m); // Symetric case ***
    //    TPZSymetricSpStructMatrix matskl(cmesh_m); // Non symetric case ***
    
    matskl.SetNumThreads(numthreads);
    analysis.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    analysis.Assemble(); // Assembla the global matrix
    
    
#ifdef PZDEBUG
//    std::ofstream filestiff("stiffness.txt");
//    an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
//    
//    std::ofstream filerhs("rhs.txt");
//    an.Rhs().Print("R = ",filerhs,EMathematicaInput);
#endif
    
    std::cout << "Solving Matrix " << std::endl;
    analysis.Solve();
    
#ifdef PZDEBUG
//    TPZFMatrix<STATE> solucao=cmesh_m->Solution();
//    std::ofstream solout("sol.txt");
//    solucao.Print("Sol",solout,EMathematicaInput);
//    
//    std::ofstream fileAlpha("alpha.txt");
//    an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
#endif
    
    // Calculation of error
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,3> Errors;
    ofstream ErroOut("Erro.txt");
    analysis.SetExact(Sol_exact);
    bool store_errors = false;
    analysis.PostProcessError(Errors, store_errors, ErroOut);
    
    
    
    // Post Process (Paraview)
    std::cout << " Post Processing " << std::endl;
    std::string plotfile("HdivDarcy.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    scalnames.Push("P_exact");
    scalnames.Push("f");
    vecnames.Push("V");
    vecnames.Push("V_exact");
    //        vecnames.Push("V_exactBC");
    
    
    int postProcessResolution = 4; //  Keep low as possible
    
    int dim = gmesh->Dimension();
    analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    analysis.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;

   
    return 0;
}



// Creating geometric mesh
TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy)
{
    
    int64_t id, index;
    int dim = 2;
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(dim);
    
    TPZVec <REAL> coord (3,0.);
    
    // Initialization of nodes:
    for(int i = 0; i < ny; i++)
    {
        for(int j = 0; j < nx; j++)
        {
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = (i)*hy/(ny - 1);
            coord[2] = 0.;
            //Get the index of the node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    // Point 1
    TPZVec<int64_t> pointtopology(1);
    pointtopology[0] = 0;
    
    gmesh->CreateGeoElement(EPoint,pointtopology,m_matPoint,id);
    
    
    // Auxiliary vector to store the connections between elements:
    TPZVec <int64_t> connect(4,0);
    
    // Connectivity of elements:
    for(int i = 0; i < (ny - 1); i++)
    {
        for(int j = 0; j < (nx - 1); j++)
        {
            index      = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,m_matID,id);
        }
    }
    
    
    // Generating neighborhood information:
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    
    int64_t el, numelements = gmesh->NElements();
    TPZManVector <int64_t> TopolPlate(4);
    
    for (el=0; el<numelements; el++)
    {
        int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
        TPZGeoEl *plate = gmesh->ElementVec()[el];
        for (int i=0; i<4; i++)
        {
            TopolPlate[i] = plate->NodeIndex(i);
        }
        
        // The boundary conditions:
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        
        // In face x = 1
        TPZVec<int64_t> ncoordzbottVec(0);
        TPZVec<int64_t> ncoordztopVec(0);
        TPZVec<int64_t> ncoordzleftVec(0);
        TPZVec<int64_t> ncoordzrightVec(0);
        
        int64_t sizeOfbottVec = 0;
        int64_t sizeOftopVec = 0;
        int64_t sizeOfleftVec = 0;
        int64_t sizeOfrightVec = 0;
        
        for (int64_t i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[1] == hy)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.)
            {
                sizeOfleftVec++;
                ncoordzleftVec.Resize(sizeOfleftVec);
                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx)
            {
                sizeOfrightVec++;
                ncoordzrightVec.Resize(sizeOfrightVec);
                ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2)
        {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,m_matBCbott);
            TPZGeoElBC(platesidebott,m_matIntBCbott);
        }
        
        if (sizeOftopVec == 2)
        {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,m_matBCtop);
            TPZGeoElBC(platesidetop,m_matIntBCtop);
        }
        
        if (sizeOfleftVec == 2)
        {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,m_matBCleft);
            TPZGeoElBC(platesideleft,m_matIntBCleft);
        }
        
        if (sizeOfrightVec == 2)
        {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,m_matBCright);
            TPZGeoElBC(platesideright,m_matIntBCright);
        }
        
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;
        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftVec.Resize(0);
        sizeOfleftVec = 0;
        ncoordzrightVec.Resize(0);
        sizeOfrightVec = 0;
        
    }
    

    
    
    // Creating interface (Generalized):
    
    TPZVec<int64_t> nodint(2);
    for(int i = 0; i < (ny - 1); i++)
    {
        for(int j = 0; j < (nx - 1); j++)
        {
            if(j>0 && j<(nx-1))
            {
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*(i+1);
                gmesh->CreateGeoElement(EOned, nodint, m_matInterface, index); // Creating interface element (GeoElement)
                
            }
            if(i>0 &&j <(ny-1))
            {
                nodint[0]=j+ny*i;
                nodint[1]=j+ny*i+1;
                gmesh->CreateGeoElement(EOned, nodint, m_matInterface, index); // Creating interface element (GeoElement)
                
            }
            
        }
    }
    
    
    id++;
    
    gmesh->AddInterfaceMaterial(m_quadmat1, m_quadmat2, m_quadmat3);
    gmesh->AddInterfaceMaterial(m_quadmat2, m_quadmat1, m_quadmat3);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    
    return gmesh;
}


TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}



TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{

    // Create computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); // Insert order of approximation polynomial
    cmesh->SetDimModel(m_dim); // Insert dimension of problem
    
    
    // Definition of approximation space
    
    TPZMat2dLin *material = new TPZMat2dLin(m_matID); //Creating material
    
    cmesh->InsertMaterialObject(material); //Inserts material into the mesh
    
    cmesh->SetAllCreateFunctionsHDiv(); //Creating HDIV Functions
        
        // Dimension of material (for HDiv):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    
    // Condition of contours:
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, m_matBCbott, m_dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond0); //Insert material no mesh
    
    TPZMaterial * BCond1 = material->CreateBC(material, m_matBCtop, m_dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond1);
    
    TPZMaterial * BCond2 = material->CreateBC(material, m_matBCleft, m_dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond2);
    
    TPZMaterial * BCond3 = material->CreateBC(material, m_matBCright, m_dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond3);
    
    
    // Creating computational elements that manage the approximation space of the mesh:
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
    }
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{

   
    // Create computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(m_dim);
    
    cmesh->SetAllCreateFunctionsContinuous(); // Create approximation H1 space
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    TPZMat2dLin *material = new TPZMat2dLin(m_matID); // Create the material
    cmesh->InsertMaterialObject(material); //Insert material in mesh
    
    // Dimension of materials (for H1 space):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    // Condition of contours
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
//        val2(0,0) = 0.0; // px -> 0
//        val2(1,0) = 0.0; // py -> 0
//    
//        TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond0);
//    
//        TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond1);
//    
//        TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond2);
//    
//        TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2);
//        cmesh->InsertMaterialObject(BCond3);
    
    
    //    Pressure point:
    //
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    TPZMaterial * BCPoint = material->CreateBC(material, m_matPoint, m_pointtype, val3, val4);
    cmesh->InsertMaterialObject(BCPoint);
    
    
    
    // Creating computational elements that manage the space of the mesh
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
    }
    std::set<int> materialids;
    materialids.insert(m_matID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    return cmesh;
    
}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    // Creating computational mesh:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(m_dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    // Create the material:
    TPZDarcyMaterial *material = new TPZDarcyMaterial(m_matID,m_dim,fviscosity,fpermeability,ftheta);
    
    // Insert material in mesh
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact);
    
    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(solp);
    cmesh->InsertMaterialObject(material);
    
    // Condition of contours
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, m_matBCbott, m_neumann, val1, val2);
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0);
    
    TPZMaterial * BCond1 = material->CreateBC(material, m_matBCtop, m_neumann, val1, val2);
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1);
    
    TPZMaterial * BCond2 = material->CreateBC(material, m_matBCleft, m_neumann, val1, val2);
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2);
    
    TPZMaterial * BCond3 = material->CreateBC(material, m_matBCright, m_neumann, val1, val2);
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3);
    
    // Points
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    val4(0,0)=0.0;
    TPZMaterial * BCPoint = material->CreateBC(material, m_matPoint, m_pointtype, val3, val4);
    cmesh->InsertMaterialObject(BCPoint);
    
    
#ifdef PZDEBUG
//        int ncel = cmesh->NElements();
//        for(int i =0; i<ncel; i++)
//        {
//            TPZCompEl * compEl = cmesh->ElementVec()[i];
//            if(!compEl) continue;
//            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//            if(facel)DebugStop();
//        }
#endif
    
    
     // Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    const REAL Pi=M_PI;
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE v_x = -2.*Pi*cos(2.*Pi*xv)*sin(2.*Pi*yv);
    STATE v_y = -2.*Pi*cos(2.*Pi*yv)*sin(2.*Pi*xv);
    STATE pressure= sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    sol[0]=v_x;
    sol[1]=v_y;
    sol[2]=pressure;
    
    // vx direction
    dsol(0,0)= 4.*Pi*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    dsol(0,1)= -4.*Pi*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    
    // vy direction
    dsol(1,0)= -4.*Pi*Pi*cos(2.*Pi*xv)*cos(2.*Pi*yv);
    dsol(1,1)= 4.*Pi*Pi*sin(2.*Pi*xv)*sin(2.*Pi*yv);
    
    // Gradient of pressure
    dsol(2,0)= -v_x;
    dsol(2,1)= -v_y;
    
    
}

void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(1);
    const REAL Pi=M_PI;
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE f_x = -8.0*Pi*Pi*sin(2.0*Pi*xv)*sin(2.0*Pi*yv);
    
    f[0] = f_x;
    
    
}



void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        if (celstack.size() != 2) {
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        int64_t index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}

