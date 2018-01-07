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
            KRATOS_CATCH("")
        }

        /**
         * Turn back information as a string.
         * (DEACTIVATED)
         */
        //std::string Info();

        /**
         * Print information about this object.
         * (DEACTIVATED)
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;

        /**
         * Print object's data.
         * (DEACTIVATED)
         */
        //virtual void PrintData(std::ostream& rOStream) const;

    protected:

        Element::Pointer pSlave() const {return mpSlaveElement;}
        Element::Pointer pMaster() const {return mpMasterElement;}

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

