/*
 *  TPZStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMatWithMem.h"
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"


#ifndef TPZSTOKESMATERIAL
#define TPZSTOKESMATERIAL



class TPZStokesMaterial : public TPZMatWithMem<TPZFMatrix<STATE>, TPZDiscontinuousGalerkin >  {
    
protected:
    
    /// dimension of the material
    int fDimension;
    
    /// Aproximation Space for velocity
    int fSpace;
    
    /// viscosidade
    STATE fViscosity;
    
    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
    STATE fk;
    
    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;
    
    STATE fSigma;
    
public:
    
    
    /**
     * Empty Constructor
     */
    TPZStokesMaterial();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZStokesMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZStokesMaterial(const TPZStokesMaterial &mat);
    
    /**
     * Destructor
     */
    ~TPZStokesMaterial();
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec);
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data);
    
    void SetPermeability(REAL perm) {
        fk = perm;
    }
    
    /** returns the name of the material */
    std::string Name() {
        return "TPZStokesMaterial";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 4;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);

    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** index of velocity */
    int VIndex(){ return 0; }
    
    /** index of pressure */
    int PIndex(){ return 1; }
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <class TVar>
    TVar Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    STATE InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);
    
    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<STATE> &GradU );
    
    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<STATE> &GradU );
    
    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    

    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    
    // Contribute Methods being used
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
   
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
   
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }

    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);
    
    
    virtual int NEvalErrors() {return 6;}
    
    /**
     * It computes errors contribution in differents spaces.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);
    
};

#endif
