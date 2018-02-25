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
//  Date:            1 Feb 2018
//


#if !defined(KRATOS_CUBIC_SPLINE_FUNCTION_H_INCLUDED )
#define  KRATOS_CUBIC_SPLINE_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/zero_function.h"
#include "custom_external_libraries/spline.h"

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
///@name  CubicSplineFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
* Implementation of cubic spline function
* REF: http://kluge.in-chemnitz.de/opensource/spline/
*/
template<int TDerivDegree>
class CubicSplineFunction : public FunctionR1R1
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CubicSplineFunction
    KRATOS_CLASS_POINTER_DEFINITION(CubicSplineFunction);

    typedef FunctionR1R1 BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CubicSplineFunction() : BaseType()
    {}

    /// Copy constructor.
    CubicSplineFunction(CubicSplineFunction const& rOther)
    : BaseType(rOther), mS(rOther.mS)
    {}

    /// Destructor.
    virtual ~CubicSplineFunction()
    {}


    virtual BaseType::Pointer CloneFunction() const
    {
        return BaseType::Pointer(new CubicSplineFunction(*this));
    }


    void SetPoints(const std::vector<double>& t, const std::vector<double>& x)
    {
        mS.set_points(t, x);
    }

    void SetLeftBoundary(const int& type, const double& value)
    {
        mS.set_left_boundary(static_cast<tk::spline::bd_type>(type), value);
    }

    void SetRightBoundary(const int& type, const double& value)
    {
        mS.set_right_boundary(static_cast<tk::spline::bd_type>(type), value);
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        if (TDerivDegree == 0)
            return mS(P);
        else if (TDerivDegree > 0)
            return mS.deriv(TDerivDegree, P);
    }


    virtual double GetDerivative(const int& component, const InputType& P) const
    {
        return mS.deriv(TDerivDegree+1, P);
    }


    virtual double GetSecondDerivative(const int& component_1, const int& component_2, const InputType& P) const
    {
        return mS.deriv(TDerivDegree+2, P);
    }


    // virtual FunctionR1R1::Pointer GetDiffFunction(const int& component) const
    // {
    //     typename CubicSplineFunction<TDerivDegree+1>::Pointer pFunc = boost::make_shared<CubicSplineFunction<TDerivDegree+1> >();
    //     pFunc->SetPoints(mS.x(), mS.y());
    //     return pFunc;
    // }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return "S";
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
        return "Cubic Spline Function";
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

    tk::spline mS;

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
    CubicSplineFunction& operator=(CubicSplineFunction const& rOther);

    ///@}

}; // Class CubicSplineFunction



template<>
class CubicSplineFunction<3> : public CubicSplineFunction<0>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(CubicSplineFunction);

    typedef FunctionR1R1 BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return mS.deriv(3, P);
    }


    virtual double GetDerivative(const int& component, const InputType& P) const
    {
        return mS.deriv(4, P);
    }


    virtual double GetSecondDerivative(const int& component_1, const int& component_2, const InputType& P) const
    {
        return mS.deriv(5, P);
    }


    // virtual FunctionR1R1::Pointer GetDiffFunction(const int& component) const
    // {
    //     return FunctionR1R1::Pointer(new ZeroFunction<FunctionR1R1>());
    // }
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream CubicSplineFunction
template<int TDerivDegree>
inline std::istream& operator >> (std::istream& rIStream, CubicSplineFunction<TDerivDegree>& rThis)
{}

/// output stream CubicSplineFunction
template<int TDerivDegree>
inline std::ostream& operator << (std::ostream& rOStream, const CubicSplineFunction<TDerivDegree>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CUBIC_SPLINE_FUNCTION_H_INCLUDED  defined
