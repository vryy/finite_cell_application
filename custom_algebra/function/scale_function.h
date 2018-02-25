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


#if !defined(KRATOS_SCALE_FUNCTION_H_INCLUDED )
#define  KRATOS_SCALE_FUNCTION_H_INCLUDED



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
///@name  ScaleFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ScaleFunction
*/
template<class TFunction>
class ScaleFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ScaleFunction
    KRATOS_CLASS_POINTER_DEFINITION(ScaleFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScaleFunction(const double a, const typename BaseType::Pointer p_func)
    : BaseType(), ma(a), mp_func(p_func)
    {}

    /// Copy constructor.
    ScaleFunction(ScaleFunction const& rOther)
    : BaseType(rOther), ma(rOther.ma), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~ScaleFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual typename BaseType::Pointer CloneFunction() const
    {
        return typename BaseType::Pointer(new ScaleFunction(*this));
    }


    virtual OutputType GetValue(const InputType& P) const
    {
        return ma*mp_func->GetValue(P);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        std::stringstream ss;
        if(ma != 1.0)
            ss << ma << "*";
        ss << mp_func->GetFormula(Format);
        return ss.str();
    }


    virtual typename BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return typename BaseType::Pointer(new ScaleFunction(ma, mp_func->GetDiffFunction(component)));
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
        return "Scale Function of " + mp_func->Info();
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

    double ma;
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
    ScaleFunction& operator=(ScaleFunction const& rOther);

    ///@}

}; // Class ScaleFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream ScaleFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, ScaleFunction<TFunction>& rThis)
{}

/// output stream ScaleFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const ScaleFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SCALE_FUNCTION_H_INCLUDED  defined
