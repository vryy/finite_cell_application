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
//  Date:            10 Feb 2017
//


#if !defined(KRATOS_LEVEL_SET_H_INCLUDED )
#define  KRATOS_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"
#include "custom_algebra/function.h"


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
/** Abstract class for a level set in space, both for implicit level set or nodal interpolated level set
*/
class LevelSet : public Function<typename Element::GeometryType::PointType::PointType, double>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LevelSet
    KRATOS_CLASS_POINTER_DEFINITION(LevelSet);

    typedef Function<typename Element::GeometryType::PointType::PointType, double> BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LevelSet() {}

    /// Destructor.
    virtual ~LevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    virtual std::size_t WorkingSpaceDimension() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// inherit from Function
    virtual double GetValue(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    virtual double GetValue(GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    virtual Vector GetGradient(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    virtual Vector GetGradient(GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    int CutStatus(Element::Pointer& p_elem) const
    {
        return CutStatus(p_elem->GetGeometry());
    }


    /// Check if a geometry is cut by the level set
    /// 0: the cell is inside the domain bounded by level set
    /// 1: outside
    /// -1: the cell is cut by level set
    int CutStatus(GeometryType& r_geom) const
    {
        std::set<std::size_t> in_list, out_list;
        for(std::size_t v = 0; v < r_geom.size(); ++v)
        {
            double phi = this->GetValue(r_geom[v]);
            if(phi <= 0.0)
                in_list.insert(v);
            else
                out_list.insert(v);
        }

        int stat;
        if(in_list.size() == r_geom.size())
        {
            stat = 0;
        }
        else
        {
            if(out_list.size() == r_geom.size())
                stat = 1;
            else
                stat = -1;
        }

        return stat;
    }


    /// Compute the intersection of the level set with a line connect by 2 points. Note that, the checking of the level set with the line is not performed. Hence one should ensure that before calling this function.
    PointType Bisect(const PointType& P1, const PointType& P2, const double& tol) const
    {
        double f1 = this->GetValue(P1);
        double f2 = this->GetValue(P2);
        if(f1*f2 > 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "Bisect does not work with two end at the same side", "")

        double left = 0.0; 
        double right = 1.0;

        bool converged = false;
        PointType P;
        while(!converged)
        {
            double mid = (left+right)/2;
            P = P1 + mid*(P2-P1);
            double fm = this->GetValue(P);

            if(fabs(fm) < tol)
            {
                converged = true;
            }
            else
            {
                if(fm*f1 < 0.0)
                {
                    right = mid;
                    f2 = fm;
                }
                else
                {
                    left = mid;
                    f1 = fm;
                }

                if(right-left < tol)
                    converged = true;
            }
        }

        return P;
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
        return "Level Set";
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
    LevelSet& operator=(LevelSet const& rOther);

    /// Copy constructor.
    LevelSet(LevelSet const& rOther);


    ///@}

}; // Class LevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, LevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const LevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LEVEL_SET_H_INCLUDED  defined
