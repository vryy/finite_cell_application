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
//  Date:            12 Apr 2017
//


#if !defined(KRATOS_MATHPRESSO_FUNCTION_H_INCLUDED )
#define  KRATOS_MATHPRESSO_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"
#include "mathpresso/mathpresso.h"


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
///@name  MathPressoFunctionR3R1
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general MathPressoFunctionR3R1
*/
class MathPressoFunctionR3R1 : public FunctionR3R1
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MathPressoFunctionR3R1
    KRATOS_CLASS_POINTER_DEFINITION(MathPressoFunctionR3R1);

    typedef FunctionR3R1 BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MathPressoFunctionR3R1(const std::string& expression)
    : mStr(expression)
    {
        this->Initialize(expression);
    }

    /// Copy constructor.
    MathPressoFunctionR3R1(MathPressoFunctionR3R1 const& rOther)
    : BaseType(rOther), mStr(rOther.mStr)
    {
        this->Initialize(rOther.mStr);
    }

    /// Destructor.
    virtual ~MathPressoFunctionR3R1()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    void Initialize(const std::string& expression)
    {
        // Initialize the context by adding MathPresso built-ins. Without this line
        // functions like round(), sin(), etc won't be available.
        mCtx.addBuiltIns();

        // Let the context know the name of the variables we will refer to and
        // their positions in the data pointer. We will use an array of 3 doubles,
        // so index them by using `sizeof(double)`, like a normal C array.
        //
        // The `addVariable()` also contains a third parameter that describes
        // variable flags, use `kVariableRO` to make a certain variable read-only.
        mCtx.addVariable("x", 0 * sizeof(double));
        mCtx.addVariable("y", 1 * sizeof(double));
        mCtx.addVariable("z", 2 * sizeof(double));

        // Compile the expression.
        //
        // The create parameters are:
        //   1. `mathpresso::Context&` - The expression's context / environment.
        //   2. `const char* body` - The expression body.
        //   3. `unsigned int` - Options, just pass `mathpresso::kNoOptions`.
        mathpresso::Error err = mExp.compile(mCtx, expression.c_str(), mathpresso::kNoOptions);

        // Handle possible syntax or compilation error.
        if (err != mathpresso::kErrorOk)
        {
            KRATOS_THROW_ERROR(std::logic_error, "Expression Error: ", err);
        }
    }


    virtual BaseType::Pointer CloneFunction() const
    {
        return BaseType::Pointer(new MathPressoFunctionR3R1(*this));
    }


    virtual double GetValue(const InputType& P) const
    {
        double data[3];
        data[0] = P[0];
        data[1] = P[1];
        data[2] = P[2];
        return mExp.evaluate(data);
    }


    virtual std::string GetFormula(const std::string& Format) const
    {
        return mStr;
    }


//    virtual typename BaseType::Pointer GetDiffFunction(const int& component) const
//    {
//        return typename BaseType::Pointer(new ZeroFunction<BaseType>());
//    }


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
        return "MathPresso Function";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << mStr;
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

    mathpresso::Context mCtx;
    mathpresso::Expression mExp;
    std::string mStr;

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
    MathPressoFunctionR3R1& operator=(MathPressoFunctionR3R1 const& rOther);

    ///@}

}; // Class MathPressoFunctionR3R1

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream MathPressoFunctionR3R1
inline std::istream& operator >> (std::istream& rIStream, MathPressoFunctionR3R1& rThis)
{}

/// output stream MathPressoFunctionR3R1
inline std::ostream& operator << (std::ostream& rOStream, const MathPressoFunctionR3R1& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MATHPRESSO_FUNCTION_H_INCLUDED  defined
