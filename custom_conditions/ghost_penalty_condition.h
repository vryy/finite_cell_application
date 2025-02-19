//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 6 Jan 2018 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_GHOST_PENALTY_CONDITION_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_CONDITION_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/legacy_structural_app_vars.h"

namespace Kratos
{
/**
 * Abstract class to represent a condition for ghost penalty stabilization.
 */
class GhostPenaltyCondition : public Condition
{
public:
    // Counted pointer of GhostPenaltyCondition
    KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyCondition);

    /**
     * Default constructor.
     */
    GhostPenaltyCondition()
    {}

    GhostPenaltyCondition( IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition( NewId, pGeometry )
    {}

    GhostPenaltyCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition( NewId, pGeometry, pProperties )
    {}

    GhostPenaltyCondition( IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           Element::Pointer pSlaveElement,
                           Element::Pointer pMasterElement,
                           PropertiesType::Pointer pProperties )
        : Condition( NewId, pGeometry, pProperties )
        , mpSlaveElement(pSlaveElement), mpMasterElement(pMasterElement)
    {}

    /**
     * Destructor.
     */
    ~GhostPenaltyCondition() override
    {}

    virtual Condition::Pointer Create(IndexType NewId,
                                      GeometryType::Pointer pGeometry,
                                      Element::Pointer pSlaveElement,
                                      Element::Pointer pMasterElement,
                                      PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer(new GhostPenaltyCondition(NewId,
                                  pGeometry, pSlaveElement, pMasterElement, pProperties));
    }

    /**
     * Operations.
     */

    /**
     * Calculates the local system contributions for this contact element
     */
    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo ) override
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo ) override
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType matrix = Matrix();
        CalculateAll( matrix, rRightHandSideVector,
                      rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag,
                      CalculateResidualVectorFlag );
    }

    void CalculateDampingMatrix( MatrixType& rDampMatrix, const ProcessInfo& rCurrentProcessInfo ) override
    {
        rDampMatrix.resize(0, 0, false);
    }

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo ) const override
    {
        rResult.resize(0);
    }

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo ) const override
    {
        ConditionalDofList.resize(0);
    }

    void Initialize( const ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_TRY

        // integration rule
        if (this->Has( INTEGRATION_ORDER ))
        {
            if (this->GetValue(INTEGRATION_ORDER) == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if (this->GetValue(INTEGRATION_ORDER) == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if (this->GetValue(INTEGRATION_ORDER) == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if (this->GetValue(INTEGRATION_ORDER) == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if (this->GetValue(INTEGRATION_ORDER) == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_ERROR << "KinematicLinear element does not support for integration order " << this->GetValue(INTEGRATION_ORDER);
        }
        else if (GetProperties().Has( INTEGRATION_ORDER ))
        {
            if (GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
            }
            else if (GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
            else if (GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
            }
            else if (GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
            }
            else if (GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
            else
                KRATOS_ERROR << "KinematicLinear element does not support for integration order " << GetProperties()[INTEGRATION_ORDER];
        }
        else
        {
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();    // default method
        }

        KRATOS_CATCH("")
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "GhostPenaltyCondition #" << Id();
        return buffer.str();
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " faces:";
        for (std::size_t i = 0; i < GetGeometry().size(); ++i)
        {
            rOStream << " " << GetGeometry()[i].Id();
        }
        rOStream << ", slave: " << mpSlaveElement->Id()
                 << ", master: " << mpMasterElement->Id()
                 << std::endl;
    }

protected:

    Element::Pointer pSlave() const {return mpSlaveElement;}
    Element::Pointer pMaster() const {return mpMasterElement;}

    IntegrationMethod mThisIntegrationMethod;

private:

    Element::Pointer mpSlaveElement;
    Element::Pointer mpMasterElement;

    friend class Serializer;

    void save ( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
    }

    void load ( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
    }

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        rLeftHandSideMatrix.resize(0, 0, false);
        rRightHandSideVector.resize(0, false);

        KRATOS_CATCH("")
    }

}; // Class GhostPenaltyCondition

}  // namespace Kratos.

#endif // KRATOS_GHOST_PENALTY_CONDITION_H_INCLUDED defined
