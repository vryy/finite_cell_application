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


#if !defined(KRATOS_PRODUCT_FUNCTION_H_INCLUDED )
#define  KRATOS_PRODUCT_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function.h"
#include "custom_algebra/sum_function.h"
#include "custom_algebra/product_function.h"


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
///@name  ProductFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ProductFunction
*/

class ProductFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ProductFunction
    KRATOS_CLASS_POINTER_DEFINITION(ProductFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ProductFunction(const BaseType::Pointer& p_func_1, const BaseType::Pointer& p_func_2)
    : mp_func_1(p_func_1), mp_func_2(p_func_2)
    {}

    /// Destructor.
    virtual ~ProductFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return mp_func_1->GetValue(P) * mp_func_2->GetValue(P);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return mp_func_1->GetFormula(Format) + "*" + mp_func_2->GetFormula(Format);
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(
                    new SumFunction(
                        BaseType::Pointer(
                            new ProductFunction(
                                mp_func_1->GetDiffFunction(component),
                                mp_func_2
                            )
                        ),
                        BaseType::Pointer(
                            new ProductFunction(
                                mp_func_1,
                                mp_func_2->GetDiffFunction(component)
                            )
                        )
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
        return "Product Function of " + mp_func_1->Info() + " and " + mp_func_2->Info();
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

    const BaseType::Pointer mp_func_1;
    const BaseType::Pointer mp_func_2;

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
    ProductFunction& operator=(ProductFunction const& rOther);

    /// Copy constructor.
    ProductFunction(ProductFunction const& rOther);


    ///@}

}; // Class ProductFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream ProductFunction
inline std::istream& operator >> (std::istream& rIStream, ProductFunction& rThis)
{}

/// output stream ProductFunction
inline std::ostream& operator << (std::ostream& rOStream, const ProductFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRODUCT_FUNCTION_H_INCLUDED  defined
