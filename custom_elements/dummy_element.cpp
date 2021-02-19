//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 28 Mar 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/dummy_element.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
DummyElement::DummyElement()
{
}

DummyElement::DummyElement( IndexType NewId,
                              GeometryType::Pointer pGeometry)
: Element( NewId, pGeometry )
{
}

DummyElement::DummyElement( IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
: Element( NewId, pGeometry, pProperties )
{
}

/**
 * Destructor. Never to be called manually
 */
DummyElement::~DummyElement()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

Element::Pointer DummyElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer DummyElement::Create(IndexType NewId, GeometryType::Pointer pGeom,
                                        PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DummyElement(NewId, pGeom, pProperties));
}

void DummyElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
void DummyElement::CalculateRightHandSide( VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
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
void DummyElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo)
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
void DummyElement::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  const ProcessInfo& rCurrentProcessInfo,
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
void DummyElement::EquationIdVector( EquationIdVectorType& rResult,
                                     const ProcessInfo& CurrentProcessInfo) const
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
void DummyElement::GetDofList( DofsVectorType& ElementalDofList,
                               const ProcessInfo& CurrentProcessInfo) const
{
    ElementalDofList.resize(0);
}

} // Namespace Kratos

