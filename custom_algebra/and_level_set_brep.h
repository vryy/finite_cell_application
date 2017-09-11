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
//  Date:            9 Sep 2017
//


#if !defined(KRATOS_AND_LEVEL_SET_ANDLevelSetBRep_H_INCLUDED )
#define  KRATOS_AND_LEVEL_SET_ANDLevelSetBRep_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"
#include "level_set/level_set.h"


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
/** BRep representing by AND of two level sets
*/
class ANDLevelSetBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ANDLevelSetBRep
    KRATOS_CLASS_POINTER_DEFINITION(ANDLevelSetBRep);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ANDLevelSetBRep(LevelSet::Pointer pLS1, LevelSet::Pointer pLS2)
    : mpLS1(pLS1), mpLS2(pLS2)
    {}

    /// Destructor.
    virtual ~ANDLevelSetBRep() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual std::size_t WorkingSpaceDimension() const
    {
        return mpLS1->WorkingSpaceDimension();
    }


    virtual std::size_t LocalSpaceDimension() const
    {
        return mpLS1->WorkingSpaceDimension();
    }

    virtual bool IsInside(const PointType& P) const
    {
        return (mpLS1->IsInside(P) && mpLS2->IsInside(P));
    }

    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(GeometryType& r_geom) const
    {
        if(mpLS1->CutStatus(r_geom) == _OUT || mpLS2->CutStatus(r_geom) == _OUT)
        {
            return _OUT;
        }
        else
        {
            if(mpLS1->CutStatus(r_geom) == _IN && mpLS2->CutStatus(r_geom) == _IN)
            {
                return _IN;
            }
            else
            {
                return _CUT;
            }
        }
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
        return "ANDLevelSetBRep";
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

    LevelSet::Pointer mpLS1;
    LevelSet::Pointer mpLS2;

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
    ANDLevelSetBRep& operator=(ANDLevelSetBRep const& rOther);

    /// Copy constructor.
    ANDLevelSetBRep(ANDLevelSetBRep const& rOther);


    ///@}

}; // Class ANDLevelSetBRep

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ANDLevelSetBRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ANDLevelSetBRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_AND_LEVEL_SET_BREP_H_INCLUDED  defined
