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
//  Date:            28 Mar 2017
//


#if !defined(KRATOS_PARAMETRIC_VOLUME_H_INCLUDED )
#define  KRATOS_PARAMETRIC_VOLUME_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
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
/** Abstract class for a parametric volume in 3D
*/
class ParametricVolume : public FunctionR3R3
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParametricVolume
    KRATOS_CLASS_POINTER_DEFINITION(ParametricVolume);

    typedef FunctionR3R3 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParametricVolume(const FunctionR3R1::Pointer p_func_x,
        const FunctionR3R1::Pointer p_func_y, const FunctionR3R1::Pointer p_func_z)
    : BaseType(), mp_func_x(p_func_x), mp_func_y(p_func_y), mp_func_z(p_func_z)
    {}

    /// Copy constructor.
    ParametricVolume(ParametricVolume const& rOther)
    : BaseType(rOther)
    , mp_func_x(rOther.mp_func_x->CloneFunction())
    , mp_func_y(rOther.mp_func_y->CloneFunction())
    , mp_func_z(rOther.mp_func_z->CloneFunction())
    {}

    /// Destructor.
    virtual ~ParametricVolume() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// inherit from Function
    virtual BaseType::Pointer CloneFunction() const
    {
        return BaseType::Pointer(new ParametricVolume(*this));
    }


    /// inherit from Function
    virtual OutputType GetValue(const InputType& T) const
    {
        OutputType P;

        P[0] = mp_func_x->GetValue(T);
        P[1] = mp_func_y->GetValue(T);
        P[2] = mp_func_z->GetValue(T);

        return P;
    }


    /// inherit from Function
    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(
                    new ParametricVolume(
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
        return "Parametric Volume";
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

    FunctionR3R1::Pointer mp_func_x;
    FunctionR3R1::Pointer mp_func_y;
    FunctionR3R1::Pointer mp_func_z;

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
    ParametricVolume& operator=(ParametricVolume const& rOther);

    ///@}

}; // Class ParametricVolume

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ParametricVolume& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ParametricVolume& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PARAMETRIC_VOLUME_H_INCLUDED  defined
