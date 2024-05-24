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

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) final;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                               VectorType& rRightHandSideVector,
                               const ProcessInfo& rCurrentProcessInfo) final;

    void CalculateRightHandSide( VectorType& rRightHandSideVector,
                                 const ProcessInfo& rCurrentProcessInfo) final;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) final;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) final;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo) const final;

    void GetDofList( DofsVectorType& ConditionalDofList,
                     const ProcessInfo& CurrentProcessInfo) const final;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

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

    void save ( Serializer& rSerializer ) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Condition )
    }

    void load ( Serializer& rSerializer ) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Condition )
    }

}; // Class ElementWrapperCondition

}  // namespace Kratos.


#endif // KRATOS_ELEMENT_WRAPPER_CONDITION_H_INCLUDED defined

