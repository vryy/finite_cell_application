//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         finite_cell_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            22 Feb 2017
//


#if !defined(KRATOS_NEGATE_FUNCTION_H_INCLUDED )
#define  KRATOS_NEGATE_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function.h"


namespace Kratos
{
///@addtogroup FiniteCellApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  NegateFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general NegateFunction
*/

class NegateFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NegateFunction
    KRATOS_CLASS_POINTER_DEFINITION(NegateFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NegateFunction(const BaseType::Pointer& p_func)
    : mp_func(p_func)
    {}

    /// Destructor.
    virtual ~NegateFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return -mp_func->GetValue(P);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        std::stringstream ss;
        ss << "-" << mp_func->GetFormula(Format);
        return ss.str();
    }


    virtual Function::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new NegateFunction(mp_func->GetDiffFunction(component)));
    }


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
    virtual std::string Info() const
    {
        return "Negate Function of " + mp_func->Info();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    const BaseType::Pointer mp_func;

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
    NegateFunction& operator=(NegateFunction const& rOther);

    /// Copy constructor.
    NegateFunction(NegateFunction const& rOther);


    ///@}

}; // Class NegateFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream NegateFunction
inline std::istream& operator >> (std::istream& rIStream, NegateFunction& rThis)
{}

/// output stream NegateFunction
inline std::ostream& operator << (std::ostream& rOStream, const NegateFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEGATE_FUNCTION_H_INCLUDED  defined
