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
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

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
        virtual ~GhostPenaltyCondition()
        {}

        virtual Condition::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeometry,
            Element::Pointer pSlaveElement,
            Element::Pointer pMasterElement,
            PropertiesType::Pointer pProperties)
        {
            return GhostPenaltyCondition::Pointer(new GhostPenaltyCondition(NewId,
                pGeometry, pSlaveElement, pMasterElement, pProperties));
        }

        /**
         * Operations.
         */

        /**
         * Calculates the local system contributions for this contact element
         */
        virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo)
        {
            //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;
            CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                          CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
        }

        virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
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

        virtual void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo)
        {
            rResult.resize(0);
        }

        virtual void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo)
        {
            ConditionalDofList.resize(0);
        }

        virtual void Initialize()
        {
            KRATOS_TRY

            Variable<int>& INTEGRATION_ORDER_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("INTEGRATION_ORDER"));

            // integration rule
            if(this->Has( INTEGRATION_ORDER_var ))
            {
                if(this->GetValue(INTEGRATION_ORDER_var) == 1)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
                }
                else if(this->GetValue(INTEGRATION_ORDER_var) == 2)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
                }
                else if(this->GetValue(INTEGRATION_ORDER_var) == 3)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
                }
                else if(this->GetValue(INTEGRATION_ORDER_var) == 4)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
                }
                else if(this->GetValue(INTEGRATION_ORDER_var) == 5)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration rule", this->GetValue(INTEGRATION_ORDER_var))
            }
            else if(GetProperties().Has( INTEGRATION_ORDER_var ))
            {
                if(GetProperties()[INTEGRATION_ORDER_var] == 1)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
                }
                else if(GetProperties()[INTEGRATION_ORDER_var] == 2)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
                }
                else if(GetProperties()[INTEGRATION_ORDER_var] == 3)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
                }
                else if(GetProperties()[INTEGRATION_ORDER_var] == 4)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
                }
                else if(GetProperties()[INTEGRATION_ORDER_var] == 5)
                {
                    mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "KinematicLinear element does not support for integration points", GetProperties()[INTEGRATION_ORDER_var])
            }
            else
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

            KRATOS_CATCH("")
        }

        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "GhostPenaltyCondition #" << Id();
            return buffer.str();
        }

        /**
         * Print information about this object.
         */
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << Info();
        }

        /**
         * Print object's data.
         */
        virtual void PrintData(std::ostream& rOStream) const
        {
            rOStream << " faces:";
            for (std::size_t i = 0; i < GetGeometry().size(); ++i)
                rOStream << " " << GetGeometry()[i].Id();
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

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
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

}; // Class GhostPenaltyCondition

}  // namespace Kratos.


#endif // KRATOS_GHOST_PENALTY_CONDITION_H_INCLUDED defined

