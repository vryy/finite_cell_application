//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 17 Mar 17 $
//   Revision:            $Revision: 1.0 $
//
//
#if !defined(KRATOS_ELEMENT_WRAPPER_CONDITION_H_INCLUDED )
#define  KRATOS_ELEMENT_WRAPPER_CONDITION_H_INCLUDED


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
 */
class ElementWrapperCondition : public Condition
{
    public:
        // Counted pointer of ElementWrapperCondition
        KRATOS_CLASS_POINTER_DEFINITION(ElementWrapperCondition);

        /**
         * Default constructor.
         */
        ElementWrapperCondition();
        ElementWrapperCondition( IndexType NewId, Element::Pointer pElement);
        ElementWrapperCondition( IndexType NewId, Element::Pointer pElement, PropertiesType::Pointer pProperties);

        /**
         * Destructor.
         */
        virtual ~ElementWrapperCondition();

        /**
         * Operations.
         */

        virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

        virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                   VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                     ProcessInfo& rCurrentProcessInfo);

        virtual void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

        virtual void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

        virtual void EquationIdVector( EquationIdVectorType& rResult,
                               ProcessInfo& rCurrentProcessInfo);

        virtual void GetDofList( DofsVectorType& ConditionalDofList,
                         ProcessInfo& CurrentProcessInfo);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

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

    private:

        Element::Pointer mpElement;

        friend class Serializer;

        virtual void save ( Serializer& rSerializer ) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
        }

}; // Class ElementWrapperCondition

}  // namespace Kratos.


#endif // KRATOS_ELEMENT_WRAPPER_CONDITION_H_INCLUDED defined

