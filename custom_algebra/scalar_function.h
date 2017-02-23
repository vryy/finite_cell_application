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
//  Date:            14 Feb 2017
//


#if !defined(KRATOS_SCALAR_FUNCTION_H_INCLUDED )
#define  KRATOS_SCALAR_FUNCTION_H_INCLUDED



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
///@name  ScalarFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ScalarFunction
*/

class ScalarFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ScalarFunction
    KRATOS_CLASS_POINTER_DEFINITION(ScalarFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScalarFunction(const double& S)
    : mS(S)
    {}

    /// Destructor.
    virtual ~ScalarFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return mS;
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        std::stringstream ss;
        ss << mS;
        return ss.str();
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ZeroFunction());
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
        return "Scalar Function";
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

    const double mS;

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
    ScalarFunction& operator=(ScalarFunction const& rOther);

    /// Copy constructor.
    ScalarFunction(ScalarFunction const& rOther);


    ///@}

}; // Class ScalarFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream ScalarFunction
inline std::istream& operator >> (std::istream& rIStream, ScalarFunction& rThis)
{}

/// output stream ScalarFunction
inline std::ostream& operator << (std::ostream& rOStream, const ScalarFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SCALAR_FUNCTION_H_INCLUDED  defined
