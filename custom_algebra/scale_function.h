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
///@name  ScaleFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ScaleFunction
*/

class ScaleFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ScaleFunction
    KRATOS_CLASS_POINTER_DEFINITION(ScaleFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScaleFunction(const double& a, const BaseType::Pointer& p_func)
    : ma(a), mp_func(p_func)
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


    virtual double GetValue(const InputType& P) const
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


    virtual Function::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ScaleFunction(ma, mp_func->GetDiffFunction(component)));
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
    ScaleFunction& operator=(ScaleFunction const& rOther);

    /// Copy constructor.
    ScaleFunction(ScaleFunction const& rOther);


    ///@}

}; // Class ScaleFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream ScaleFunction
inline std::istream& operator >> (std::istream& rIStream, ScaleFunction& rThis)
{}

/// output stream ScaleFunction
inline std::ostream& operator << (std::ostream& rOStream, const ScaleFunction& rThis)
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
