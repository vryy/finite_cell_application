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
#include "custom_algebra/function/function.h"


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
template<class TFunction>
class NegateFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NegateFunction
    KRATOS_CLASS_POINTER_DEFINITION(NegateFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NegateFunction(const typename BaseType::Pointer p_func)
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


    virtual typename BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return typename BaseType::Pointer(new NegateFunction(mp_func->GetDiffFunction(component)));
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

    const typename BaseType::Pointer mp_func;

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
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, NegateFunction<TFunction>& rThis)
{}

/// output stream NegateFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const NegateFunction<TFunction>& rThis)
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
