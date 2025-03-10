//  see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Feb 7, 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED)
#define KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_components.h"
#include "includes/kratos_application.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

/** Application for finite cell and CutFEM discretization method
*/
class KRATOS_API(FINITE_CELL_APPLICATION) KratosFiniteCellApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMultiphaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFiniteCellApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFiniteCellApplication();

    /// Destructor.
    virtual ~KratosFiniteCellApplication() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Base application for Finite-Cell-based analysis";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "in KratosFiniteCellApplication:";
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    /// Assignment operator.
    KratosFiniteCellApplication& operator=(KratosFiniteCellApplication const& rOther);

    /// Copy constructor.
    KratosFiniteCellApplication(KratosFiniteCellApplication const& rOther);


    ///@}

}; // Class KratosFiniteCellApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


} // namespace Kratos

#endif // KRATOS_FINITE_CELL_APPLICATION_H_INCLUDED defined

