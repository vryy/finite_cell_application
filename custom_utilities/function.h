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
//  Date:            13 Feb 2017
//


#if !defined(KRATOS_FUNCTION_H_INCLUDED )
#define  KRATOS_FUNCTION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


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
/** Abstract class for a general function R^m->R^n
*/
template<typename TInputType, typename TOutputType>
class Function
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Function
    KRATOS_CLASS_POINTER_DEFINITION(Function);

    typedef TInputType InputType;

    typedef TOutputType OutputType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Function() {}

    /// Destructor.
    virtual ~Function() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual TOutputType GetValue(const TInputType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Call the base class", __FUNCTION__)
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
        return "Function R^m->R^n";
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
    Function& operator=(Function const& rOther);

    /// Copy constructor.
    Function(Function const& rOther);


    ///@}

}; // Class Function

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<typename TInputType, typename TOutputType>
inline std::istream& operator >> (std::istream& rIStream,
                Function<TInputType, TOutputType>& rThis)
{}

/// output stream function
template<typename TInputType, typename TOutputType>
inline std::ostream& operator << (std::ostream& rOStream,
                const Function<TInputType, TOutputType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FUNCTION_H_INCLUDED  defined
