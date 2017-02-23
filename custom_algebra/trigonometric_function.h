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
//  Date:            15 Feb 2017
//


#if !defined(KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED )
#define  KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function.h"
#include "custom_algebra/product_function.h"
#include "custom_algebra/negate_function.h"
#include "custom_algebra/pow_function.h"
#include "custom_algebra/scalar_function.h"

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
///@name  SinFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general SinFunction
*/

class CosFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CosFunction
    KRATOS_CLASS_POINTER_DEFINITION(CosFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CosFunction(const BaseType::Pointer& p_func)
    : mp_func(p_func)
    {}

    /// Destructor.
    virtual ~CosFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return cos(mp_func->GetValue(P));
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return "cos(" + mp_func->GetFormula(Format) + ")";
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const;

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
        return "Cos Function";
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
    CosFunction& operator=(CosFunction const& rOther);

    /// Copy constructor.
    CosFunction(CosFunction const& rOther);


    ///@}

}; // Class CosFunction

class SinFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SinFunction
    KRATOS_CLASS_POINTER_DEFINITION(SinFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SinFunction(const BaseType::Pointer& p_func)
    : mp_func(p_func)
    {}

    /// Destructor.
    virtual ~SinFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return sin(mp_func->GetValue(P));
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return "sin(" + mp_func->GetFormula(Format) + ")";
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const;


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
        return "Sin Function";
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
    SinFunction& operator=(SinFunction const& rOther);

    /// Copy constructor.
    SinFunction(SinFunction const& rOther);


    ///@}

}; // Class SinFunction

class AcosFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AcosFunction
    KRATOS_CLASS_POINTER_DEFINITION(AcosFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AcosFunction(const BaseType::Pointer& p_func)
    : mp_func(p_func)
    {}

    /// Destructor.
    virtual ~AcosFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return acos(mp_func->GetValue(P));
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return "acos(" + mp_func->GetFormula(Format) + ")";
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const;


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
        return "Acos Function";
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
    AcosFunction& operator=(AcosFunction const& rOther);

    /// Copy constructor.
    AcosFunction(AcosFunction const& rOther);


    ///@}

}; // Class AcosFunction

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream SinFunction
inline std::istream& operator >> (std::istream& rIStream, SinFunction& rThis)
{}

/// output stream SinFunction
inline std::ostream& operator << (std::ostream& rOStream, const SinFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/// input stream CosFunction
inline std::istream& operator >> (std::istream& rIStream, CosFunction& rThis)
{}

/// output stream CosFunction
inline std::ostream& operator << (std::ostream& rOStream, const CosFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/// input stream AcosFunction
inline std::istream& operator >> (std::istream& rIStream, AcosFunction& rThis)
{}

/// output stream CosFunction
inline std::ostream& operator << (std::ostream& rOStream, const AcosFunction& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block


CosFunction::BaseType::Pointer CosFunction::GetDiffFunction(const int& component) const
{
    return BaseType::Pointer(
                new NegateFunction(
                    BaseType::Pointer(
                        new ProductFunction(
                            BaseType::Pointer(new SinFunction(mp_func)),
                            mp_func->GetDiffFunction(component)
                        )
                    )
                )
            );
}

SinFunction::BaseType::Pointer SinFunction::GetDiffFunction(const int& component) const
{
    return BaseType::Pointer(
                new ProductFunction(
                    BaseType::Pointer(new CosFunction(mp_func)),
                    mp_func->GetDiffFunction(component)
                )
            );
}

AcosFunction::BaseType::Pointer AcosFunction::GetDiffFunction(const int& component) const
{
    return BaseType::Pointer(
                new NegateFunction(
                    BaseType::Pointer(
                        new ProductFunction(
                            mp_func->GetDiffFunction(component),
                            BaseType::Pointer(
                                new PowFunction(
                                    -0.5,
                                    BaseType::Pointer(
                                        new SumFunction(
                                            BaseType::Pointer( new ScalarFunction(1.0) ),
                                            BaseType::Pointer(
                                                new NegateFunction(
                                                    BaseType::Pointer( new PowFunction(2.0, mp_func) )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            );
}

}  // namespace Kratos.

#endif // KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED  defined
