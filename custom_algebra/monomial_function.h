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


#if !defined(KRATOS_MONOMIAL_FUNCTION_H_INCLUDED )
#define  KRATOS_MONOMIAL_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function.h"
#include "custom_algebra/scale_function.h"
#include "custom_algebra/zero_function.h"
#include "custom_algebra/level_set.h"


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
///@name  MonomialFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general MonomialFunction
*/

template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
std::string GetFormulaMonomialFunction(const std::string& Format)
{
    if(Format == "matlab")
    {
        std::stringstream ss;
        bool trail = false;
        if(TDegreeX != 0)
        {
            if(TDegreeX == 1)
                ss << "x";
            else
                ss << "x^" << TDegreeX;
            trail = true;
        }
        if(TDegreeY != 0)
        {
            if(trail) ss << "*";
            if(TDegreeY == 1)
                ss << "y";
            else
                ss << "y^" << TDegreeY;
            trail = true;
        }
        if(TDegreeZ != 0)
        {
            if(trail) ss << "*";
            if(TDegreeZ == 1)
                ss << "z";
            else
                ss << "z^" << TDegreeZ;
        }
        if((TDegreeX == 0) && (TDegreeY == 0) && (TDegreeZ == 0))
            ss << "1.0";
        return ss.str();
    }
}

template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
class MonomialFunction : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonomialFunction()
    {}

    /// Destructor.
    virtual ~MonomialFunction()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[1], TDegreeY) * pow(P[2], TDegreeZ);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<TDegreeX, TDegreeY, TDegreeZ>(Format);
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeX, BaseType::Pointer(new MonomialFunction<TDegreeX-1, TDegreeY, TDegreeZ>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeY, BaseType::Pointer(new MonomialFunction<TDegreeX, TDegreeY-1, TDegreeZ>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeY, BaseType::Pointer(new MonomialFunction<TDegreeX, TDegreeY, TDegreeZ-1>())));
        }
        else
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
        return "Monomial Function";
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
    MonomialFunction& operator=(MonomialFunction const& rOther);

    /// Copy constructor.
    MonomialFunction(MonomialFunction const& rOther);


    ///@}

}; // Class MonomialFunction


template<std::size_t TDegreeY, std::size_t TDegreeZ>
class MonomialFunction<0, TDegreeY, TDegreeZ> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[1], TDegreeY) * pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<0, TDegreeY, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeY, BaseType::Pointer(new MonomialFunction<0, TDegreeY-1, TDegreeZ>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeZ, BaseType::Pointer(new MonomialFunction<0, TDegreeY, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<std::size_t TDegreeX, std::size_t TDegreeZ>
class MonomialFunction<TDegreeX, 0, TDegreeZ> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<TDegreeX, 0, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeX, BaseType::Pointer(new MonomialFunction<TDegreeX-1, 0, TDegreeZ>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeZ, BaseType::Pointer(new MonomialFunction<TDegreeX, 0, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<std::size_t TDegreeX, std::size_t TDegreeY>
class MonomialFunction<TDegreeX, TDegreeY, 0> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[1], TDegreeY);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<TDegreeX, TDegreeY, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeX, BaseType::Pointer(new MonomialFunction<TDegreeX-1, TDegreeY, 0>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeY, BaseType::Pointer(new MonomialFunction<TDegreeX, TDegreeY-1, 0>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<std::size_t TDegreeX>
class MonomialFunction<TDegreeX, 0, 0> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<TDegreeX, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeX, BaseType::Pointer(new MonomialFunction<TDegreeX-1, 0, 0>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<std::size_t TDegreeY>
class MonomialFunction<0, TDegreeY, 0> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[1], TDegreeY);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<0, TDegreeY, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeY, BaseType::Pointer(new MonomialFunction<0, TDegreeY-1, 0>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<std::size_t TDegreeZ>
class MonomialFunction<0, 0, TDegreeZ> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<0, 0, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction(TDegreeZ, BaseType::Pointer(new MonomialFunction<0, 0, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction());
    }
};

template<>
class MonomialFunction<0, 0, 0> : public Function<Element::GeometryType::PointType::PointType, double>
{
public:
    /// Pointer definition of MonomialFunction
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunction);

    typedef Function<Element::GeometryType::PointType::PointType, double> BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual double GetValue(const InputType& P) const
    {
        return 1.0;
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunction<0, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ZeroFunction());
    }
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream MonomialFunction
template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
inline std::istream& operator >> (std::istream& rIStream,
            MonomialFunction<TDegreeX, TDegreeY, TDegreeZ>& rThis)
{}

/// output stream MonomialFunction
template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
inline std::ostream& operator << (std::ostream& rOStream,
            const MonomialFunction<TDegreeX, TDegreeY, TDegreeZ>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MONOMIAL_FUNCTION_H_INCLUDED  defined
