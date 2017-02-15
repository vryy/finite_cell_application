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
#include "custom_utilities/function.h"
#include "custom_utilities/level_set.h"


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

template<std::size_t DegreeX, std::size_t DegreeY, std::size_t DegreeZ>
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
        return pow(P[0], DegreeX) * pow(P[1], DegreeY) * pow(P[2], DegreeZ);
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

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream MonomialFunction
template<std::size_t DegreeX, std::size_t DegreeY, std::size_t DegreeZ>
inline std::istream& operator >> (std::istream& rIStream,
            MonomialFunction<DegreeX, DegreeY, DegreeZ>& rThis)
{}

/// output stream MonomialFunction
template<std::size_t DegreeX, std::size_t DegreeY, std::size_t DegreeZ>
inline std::ostream& operator << (std::ostream& rOStream,
            const MonomialFunction<DegreeX, DegreeY, DegreeZ>& rThis)
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
