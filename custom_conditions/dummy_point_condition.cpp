//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Mar 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes 

// External includes 

// Project includes 
#include "custom_conditions/dummy_point_condition.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyPointCondition::DummyPointCondition()
{
}

DummyPointCondition::DummyPointCondition( IndexType NewId, 
                              GeometryType::Pointer pGeometry)
: Condition( NewId, pGeometry )
{
}

DummyPointCondition::DummyPointCondition( IndexType NewId, 
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Condition( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyPointCondition::~DummyPointCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Condition::Pointer DummyPointCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyPointCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Condition::Pointer DummyPointCondition::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DummyPointCondition(NewId, pGeom, pProperties));
}

void DummyPointCondition::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}

//************************************************************************************ 
//************************************************************************************
//************************************************************************************
//************************************************************************************
/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void DummyPointCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType matrix = Matrix();
    CalculateAll( matrix, rRightHandSideVector, 
                  rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, 
                  CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this contact element's local contributions
 */
void DummyPointCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                          VectorType& rRightHandSideVector, 
                                          ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}
    //************************************************************************************
//************************************************************************************    /
/**
 * This function calculates all system contributions due to the contact problem
 * with regard to the current master and slave partners.
 * All Conditions are assumed to be defined in 2D/3D space and having 2/3 DOFs per node 
 */
void DummyPointCondition::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                  VectorType& rRightHandSideVector, 
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    rLeftHandSideMatrix.resize(0, 0, false);
    rRightHandSideVector.resize(0, false);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.    
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void DummyPointCondition::EquationIdVector( EquationIdVectorType& rResult, 
                                      ProcessInfo& CurrentProcessInfo)
{
    rResult.resize(0);
}

//************************************************************************************
//************************************************************************************
/**
 * Setting up the DOF list for the current partners.
 * All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per Node.
 * All DOF are given Master first, Slave second
 */
//************************************************************************************
//************************************************************************************
void DummyPointCondition::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    ConditionalDofList.resize(0);
}

} // Namespace Kratos

