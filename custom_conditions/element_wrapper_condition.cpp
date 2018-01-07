//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 17 Mar 2017$
//   Revision:            $Revision: 1.0 $
//
//
// System includes

// External includes

// Project includes
#include "custom_conditions/element_wrapper_condition.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
ElementWrapperCondition::ElementWrapperCondition()
{
}

ElementWrapperCondition::ElementWrapperCondition( IndexType NewId,
                              Element::Pointer pElement)
: Condition( NewId, pElement->pGetGeometry() ), mpElement(pElement)
{
}

ElementWrapperCondition::ElementWrapperCondition( IndexType NewId,
                              Element::Pointer pElement,
                              PropertiesType::Pointer pProperties)
: Condition( NewId, pElement->pGetGeometry(), pProperties ), mpElement(pElement)
{
}

/**
 * Destructor. Never to be called manually
 */
ElementWrapperCondition::~ElementWrapperCondition()
{
}


//********************************************************
//**** Operations ****************************************
//********************************************************

void ElementWrapperCondition::Initialize()
{
    mpElement->Initialize();
}

//************************************************************************************
//************************************************************************************

void ElementWrapperCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    mpElement->InitializeSolutionStep(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void ElementWrapperCondition::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpElement->InitializeNonLinearIteration(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void ElementWrapperCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    mpElement->CalculateRightHandSide( rRightHandSideVector, rCurrentProcessInfo );
}

//************************************************************************************
//************************************************************************************
/**
 * calculates this element's local contributions
 */
void ElementWrapperCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo)
{
//    std::cout << "---------------------start condition " << Id() << "---------------------------" << std::endl;
//    std::cout << "Condition " << Id() << " wraps element " << mpElement->Id() << ":";
//    for(std::size_t i = 0; i < mpElement->GetGeometry().size(); ++i)
//        std::cout << " " << mpElement->GetGeometry()[i].Id();
//    std::cout << std::endl;
    mpElement->CalculateLocalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );
//    KRATOS_WATCH(rRightHandSideVector)
//    KRATOS_WATCH(rLeftHandSideMatrix)
//    std::cout << "----------------------end condition " << Id() << "--------------------------" << std::endl;
}

//************************************************************************************
//************************************************************************************

void ElementWrapperCondition::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    mpElement->FinalizeNonLinearIteration(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void ElementWrapperCondition::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    mpElement->FinalizeSolutionStep(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
/**
* Setting up the EquationIdVector for the current partners.
* All conditions are assumed to be defined in 2D/3D space with 2/3 DOFs per node.
* All Equation IDs are given Master first, Slave second
*/
void ElementWrapperCondition::EquationIdVector( EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo)
{
    mpElement->EquationIdVector( rResult, rCurrentProcessInfo );
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
void ElementWrapperCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo)
{
    mpElement->GetDofList( rConditionalDofList, rCurrentProcessInfo );
}

} // Namespace Kratos

