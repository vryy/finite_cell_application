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
//  Date:            5 Jan 2017
//


#if !defined(KRATOS_UNION_LEVEL_SET_H_INCLUDED )
#define  KRATOS_UNION_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for union of two level sets
 * REF: Massing et al, CutFEM: Discretizing geometry and partial differential equations
*/
class UnionLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of UnionLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(UnionLevelSet);

    typedef LevelSet BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UnionLevelSet(const BaseType::Pointer p_level_set_1, const BaseType::Pointer p_level_set_2)
    : mp_level_set_1(p_level_set_1), mp_level_set_2(p_level_set_2)
    {}

    /// Destructor.
    virtual ~UnionLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual std::size_t WorkingSpaceDimension() const
    {
        return mp_level_set_1->WorkingSpaceDimension();
    }


    virtual double GetValue(const PointType& P) const
    {
        return std::min(mp_level_set_1->GetValue(P), mp_level_set_2->GetValue(P));
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        if (mp_level_set_1->GetValue(P) < mp_level_set_2->GetValue(P))
            return mp_level_set_1->GetGradient(P);
        else
            return mp_level_set_2->GetGradient(P);
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
        return "Union Level Set";
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


    const BaseType::Pointer mp_level_set_1;
    const BaseType::Pointer mp_level_set_2;


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
    UnionLevelSet& operator=(UnionLevelSet const& rOther);

    /// Copy constructor.
    UnionLevelSet(UnionLevelSet const& rOther);


    ///@}

}; // Class UnionLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, UnionLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const UnionLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_UNION_LEVEL_SET_H_INCLUDED  defined