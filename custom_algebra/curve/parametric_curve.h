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
//  Date:            24 Feb 2017
//


#if !defined(KRATOS_PARAMETRIC_CURVE_H_INCLUDED )
#define  KRATOS_PARAMETRIC_CURVE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"
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
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Abstract class for a parametric curve in 3D
*/
class ParametricCurve : public FunctionR1R3
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParametricCurve
    KRATOS_CLASS_POINTER_DEFINITION(ParametricCurve);

    typedef FunctionR1R3 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParametricCurve(const FunctionR1R1::Pointer& p_func_x,
        const FunctionR1R1::Pointer& p_func_y, const FunctionR1R1::Pointer& p_func_z)
    : BaseType(), mp_func_x(p_func_x), mp_func_y(p_func_y), mp_func_z(p_func_z)
    {}

    /// Destructor.
    virtual ~ParametricCurve() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// inherit from Function
    virtual OutputType GetValue(const InputType& t) const
    {
        OutputType P;

        P[0] = mp_func_x->GetValue(t);
        P[1] = mp_func_y->GetValue(t);
        P[2] = mp_func_z->GetValue(t);

        return P;
    }


    /// inherit from Function
    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(
                    new ParametricCurve(
                        mp_func_x->GetDiffFunction(component),
                        mp_func_y->GetDiffFunction(component),
                        mp_func_z->GetDiffFunction(component)
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
        return "Parametric Curve";
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

    FunctionR1R1::Pointer mp_func_x;
    FunctionR1R1::Pointer mp_func_y;
    FunctionR1R1::Pointer mp_func_z;

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
    ParametricCurve& operator=(ParametricCurve const& rOther);

    /// Copy constructor.
    ParametricCurve(ParametricCurve const& rOther);


    ///@}

}; // Class ParametricCurve

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ParametricCurve& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ParametricCurve& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXPLICIT_CURVE_H_INCLUDED  defined
