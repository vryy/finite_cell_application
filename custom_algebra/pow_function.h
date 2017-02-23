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


#if !defined(KRATOS_POW_FUNCTION_H_INCLUDED )
#define  KRATOS_POW_FUNCTION_H_INCLUDED



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
///@name  PowFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general PowFunction
*/

class PowFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowFunction
    KRATOS_CLASS_POINTER_DEFINITION(PowFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PowFunction(const double& a, const BaseType::Pointer& p_func)
    : ma(a), mp_func(p_func)
    {}

    PowFunction(const BaseType::Pointer& p_func, const double& a)
    : ma(a), mp_func(p_func)
    {}

    /// Destructor.
    virtual ~PowFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return pow(mp_func->GetValue(P), ma);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        std::stringstream ss;
        if(Format == "matlab")
        {
            if(ma == 1.0)
                ss << mp_func->GetFormula(Format);
            else if(ma == 0.0)
                ss << "1.0";
            else
                ss << "(" << mp_func->GetFormula(Format) << ")^" << ma;
        }
        return ss.str();
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(ma == 1.0)
            return mp_func->GetDiffFunction(component);
        else if(ma == 0.0)
            return BaseType::Pointer(new ZeroFunction());
        else
            return BaseType::Pointer(
                        new ProductFunction(
                            BaseType::Pointer(
                                new ScaleFunction(ma,
                                    BaseType::Pointer(new PowFunction(ma-1, mp_func))
                                )
                            ),
                            mp_func->GetDiffFunction(component)
                        )
                    );
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
    PowFunction& operator=(PowFunction const& rOther);

    /// Copy constructor.
    PowFunction(PowFunction const& rOther);


    ///@}

}; // Class PowFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream PowFunction
inline std::istream& operator >> (std::istream& rIStream, PowFunction& rThis)
{}

/// output stream PowFunction
inline std::ostream& operator << (std::ostream& rOStream, const PowFunction& rThis)
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
