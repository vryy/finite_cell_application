//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 1 Mar 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_DUMMY_CONDITION_H_INCLUDED )
#define  KRATOS_DUMMY_CONDITION_H_INCLUDED


// External includes

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 */
class DummyCondition : public Condition
{
    public:
        // Counted pointer of DummyCondition
        KRATOS_CLASS_POINTER_DEFINITION(DummyCondition);

        /**
         * Default constructor.
         */
        DummyCondition();
        DummyCondition( IndexType NewId, GeometryType::Pointer pGeometry);
        DummyCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~DummyCondition();

        /**
         * Operations.
         */

        Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const final;

        Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const final;

        /**
         * Calculates the local system contributions for this contact element
         */
        void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) final;

        void CalculateDampingMatrix( MatrixType& rDampingMatrix,
                                     const ProcessInfo& rCurrentProcessInfo ) final;

        void CalculateMassMatrix( MatrixType& rMassMatrix,
                                  const ProcessInfo& rCurrentProcessInfo ) final;

        void EquationIdVector( EquationIdVectorType& rResult,
                               const ProcessInfo& rCurrentProcessInfo ) const final;

        void GetDofList( DofsVectorType& ConditionalDofList,
                         const ProcessInfo& CurrentProcessInfo ) const final;

        void Initialize(const ProcessInfo& rCurrentProcessInfo);
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


    private:

        friend class Serializer;

        void save ( Serializer& rSerializer ) const final
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        void load ( Serializer& rSerializer ) final
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

        void CalculateAll( MatrixType& rLeftHandSideMatrix,
                           VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag);

}; // Class DummyCondition

}  // namespace Kratos.


#endif // KRATOS_DUMMY_CONDITION_H_INCLUDED defined

