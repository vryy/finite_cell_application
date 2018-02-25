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
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/scale_function.h"
#include "custom_algebra/function/zero_function.h"
#include "custom_algebra/level_set/level_set.h"


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
///@name  MonomialFunctionR3R1s
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general MonomialFunctionR3R1
*/

template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
std::string GetFormulaMonomialFunctionR3R1(const std::string& Format)
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
class MonomialFunctionR3R1 : public FunctionR3R1
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonomialFunctionR3R1()
    {}

    /// Copy constructor.
    MonomialFunctionR3R1(MonomialFunctionR3R1 const& rOther) : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~MonomialFunctionR3R1()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }


    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[1], TDegreeY) * pow(P[2], TDegreeZ);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, TDegreeY, TDegreeZ>(Format);
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX-1, TDegreeY, TDegreeZ>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX, TDegreeY-1, TDegreeZ>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX, TDegreeY, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
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
        return "Monomial Function R3R1";
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
    MonomialFunctionR3R1& operator=(MonomialFunctionR3R1 const& rOther);

    ///@}

}; // Class MonomialFunctionR3R1


template<std::size_t TDegreeY, std::size_t TDegreeZ>
class MonomialFunctionR3R1<0, TDegreeY, TDegreeZ> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[1], TDegreeY) * pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, TDegreeY, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR3R1<0, TDegreeY-1, TDegreeZ>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeZ, BaseType::Pointer(new MonomialFunctionR3R1<0, TDegreeY, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeX, std::size_t TDegreeZ>
class MonomialFunctionR3R1<TDegreeX, 0, TDegreeZ> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, 0, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX-1, 0, TDegreeZ>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeZ, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX, 0, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeX, std::size_t TDegreeY>
class MonomialFunctionR3R1<TDegreeX, TDegreeY, 0> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[1], TDegreeY);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, TDegreeY, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX-1, TDegreeY, 0>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX, TDegreeY-1, 0>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeX>
class MonomialFunctionR3R1<TDegreeX, 0, 0> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR3R1<TDegreeX-1, 0, 0>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeY>
class MonomialFunctionR3R1<0, TDegreeY, 0> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[1], TDegreeY);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, TDegreeY, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR3R1<0, TDegreeY-1, 0>())));
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeZ>
class MonomialFunctionR3R1<0, 0, TDegreeZ> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[2], TDegreeZ);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, 0, TDegreeZ>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 2)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeZ, BaseType::Pointer(new MonomialFunctionR3R1<0, 0, TDegreeZ-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<>
class MonomialFunctionR3R1<0, 0, 0> : public FunctionR3R1
{
public:
    /// Pointer definition of MonomialFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR3R1::Pointer CloneFunction() const
    {
        return FunctionR3R1::Pointer(new MonomialFunctionR3R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return 1.0;
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


template<std::size_t TDegreeX, std::size_t TDegreeY>
class MonomialFunctionR2R1 : public FunctionR2R1
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonomialFunctionR2R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR2R1);

    typedef FunctionR2R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonomialFunctionR2R1()
    {}

    /// Copy constructor.
    MonomialFunctionR2R1(MonomialFunctionR2R1 const& rOther) : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~MonomialFunctionR2R1()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual FunctionR2R1::Pointer CloneFunction() const
    {
        return FunctionR2R1::Pointer(new MonomialFunctionR2R1());
    }


    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX) * pow(P[1], TDegreeY);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, TDegreeY, 0>(Format);
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR2R1<TDegreeX-1, TDegreeY>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR2R1<TDegreeX, TDegreeY-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
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
        return "Monomial Function R2R1";
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
    MonomialFunctionR2R1& operator=(MonomialFunctionR2R1 const& rOther);

    ///@}

}; // Class MonomialFunctionR2R1

template<std::size_t TDegreeX>
class MonomialFunctionR2R1<TDegreeX, 0> : public FunctionR2R1
{
public:
    /// Pointer definition of MonomialFunctionR2R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR2R1);

    typedef FunctionR2R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR2R1::Pointer CloneFunction() const
    {
        return FunctionR2R1::Pointer(new MonomialFunctionR2R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[0], TDegreeX);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegreeX, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeX, BaseType::Pointer(new MonomialFunctionR2R1<TDegreeX-1, 0>())));
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<std::size_t TDegreeY>
class MonomialFunctionR2R1<0, TDegreeY> : public FunctionR2R1
{
public:
    /// Pointer definition of MonomialFunctionR2R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR2R1);

    typedef FunctionR2R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR2R1::Pointer CloneFunction() const
    {
        return FunctionR2R1::Pointer(new MonomialFunctionR2R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return pow(P[1], TDegreeY);
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, TDegreeY, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ZeroFunction<BaseType>());
        }
        else if(component == 1)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegreeY, BaseType::Pointer(new MonomialFunctionR2R1<0, TDegreeY-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

template<>
class MonomialFunctionR2R1<0, 0> : public FunctionR2R1
{
public:
    /// Pointer definition of MonomialFunctionR2R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR2R1);

    typedef FunctionR2R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR2R1::Pointer CloneFunction() const
    {
        return FunctionR2R1::Pointer(new MonomialFunctionR2R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return 1.0;
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


template<std::size_t TDegree>
class MonomialFunctionR1R1 : public FunctionR1R1
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonomialFunctionR1R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR1R1);

    typedef FunctionR1R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonomialFunctionR1R1()
    {}

    /// Destructor.
    virtual ~MonomialFunctionR1R1()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual FunctionR1R1::Pointer CloneFunction() const
    {
        return FunctionR1R1::Pointer(new MonomialFunctionR1R1());
    }


    virtual double GetValue(const InputType& P) const
    {
        return pow(P, TDegree);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<TDegree, 0, 0>(Format);
    }


    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if(component == 0)
        {
            return BaseType::Pointer(new ScaleFunction<BaseType>(TDegree, BaseType::Pointer(new MonomialFunctionR1R1<TDegree-1>())));
        }
        else
            return BaseType::Pointer(new ZeroFunction<BaseType>());
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
        return "Monomial Function R1R1";
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
    MonomialFunctionR1R1& operator=(MonomialFunctionR1R1 const& rOther);

    /// Copy constructor.
    MonomialFunctionR1R1(MonomialFunctionR1R1 const& rOther);


    ///@}

}; // Class MonomialFunctionR1R1

template<>
class MonomialFunctionR1R1<0> : public FunctionR1R1
{
public:
    /// Pointer definition of MonomialFunctionR1R1
    KRATOS_CLASS_POINTER_DEFINITION(MonomialFunctionR1R1);

    typedef FunctionR1R1 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    virtual FunctionR1R1::Pointer CloneFunction() const
    {
        return FunctionR1R1::Pointer(new MonomialFunctionR1R1());
    }

    virtual double GetValue(const InputType& P) const
    {
        return 1.0;
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        return GetFormulaMonomialFunctionR3R1<0, 0, 0>(Format);
    }

    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        return BaseType::Pointer(new ZeroFunction<BaseType>());
    }
};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream MonomialFunctionR3R1
template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
inline std::istream& operator >> (std::istream& rIStream,
            MonomialFunctionR3R1<TDegreeX, TDegreeY, TDegreeZ>& rThis)
{}

/// output stream MonomialFunctionR3R1
template<std::size_t TDegreeX, std::size_t TDegreeY, std::size_t TDegreeZ>
inline std::ostream& operator << (std::ostream& rOStream,
            const MonomialFunctionR3R1<TDegreeX, TDegreeY, TDegreeZ>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

/// input stream MonomialFunctionR1R1
template<std::size_t TDegree>
inline std::istream& operator >> (std::istream& rIStream,
            MonomialFunctionR1R1<TDegree>& rThis)
{}

/// output stream MonomialFunctionR1R1
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream,
            const MonomialFunctionR1R1<TDegree>& rThis)
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
