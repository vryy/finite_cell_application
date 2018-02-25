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
//  Date:            13 Mar 2017
//


#if !defined(KRATOS_BREP_H_INCLUDED )
#define  KRATOS_BREP_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"


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
/** Abstract class for a boundary representation
*/
class BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRep
    KRATOS_CLASS_POINTER_DEFINITION(BRep);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    static const int _CUT = -1;
    static const int _IN  = 0;
    static const int _OUT = 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRep() : mTOL(1.0e-10) {}

    /// Copy constructor.
    BRep(BRep const& rOther) : mTOL(rOther.mTOL) {}

    /// Destructor.
    virtual ~BRep() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Clone this BRep
    virtual BRep::Pointer CloneBRep() const
    {
        return BRep::Pointer(new BRep(*this));
    }

    /// Set for geometric tolerance
    void SetTolerance(const double& TOL) {mTOL = TOL;}

    /// Get for geometric tolerance
    const double GetTolerance() const {return mTOL;}

    /// Get working space dimension
    virtual std::size_t WorkingSpaceDimension() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Get local space dimension
    virtual std::size_t LocalSpaceDimension() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Check if a point is inside/outside of the BRep
    virtual bool IsInside(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Check if a point is on the boundary within a tolerance
    virtual bool IsOnBoundary(const PointType& P, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Check if an element is cut by the brep
    int CutStatus(Element::Pointer p_elem) const
    {
        return this->CutStatus(p_elem->GetGeometry());
    }


    /// Check if a geometry is cut by the brep
    int CutStatus(GeometryType::Pointer p_geom) const
    {
        return this->CutStatus(*p_geom);
    }


    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(GeometryType& r_geom) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Check if a set of points is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(const std::vector<PointType>& r_points) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }


    /// Check if an element is cut by the brep by sampling the elemental geometry
    int CutStatusBySampling(Element::Pointer p_elem, const std::size_t& nsampling) const
    {
        return this->CutStatusBySampling(p_elem->GetGeometry(), nsampling);
    }


    /// Check if a geometry is cut by the brep by sampling the geometry
    int CutStatusBySampling(GeometryType::Pointer p_geom, const std::size_t& nsampling) const
    {
        return this->CutStatusBySampling(*p_geom, nsampling);
    }


    /// Check if a geometry is cut by the level set by sampling the geometry
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatusBySampling(GeometryType& r_geom, const std::size_t& nsampling) const
    {
        std::vector<PointType> SamplingPoints;
        CreateSamplingPoints(SamplingPoints, r_geom, nsampling);
        return this->CutStatus(SamplingPoints);
    }


    /// Compute the intersection of the BRep with a line connect by 2 points.
    virtual PointType Bisect(const PointType& P1, const PointType& P2, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
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
        return "BRep";
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

    double mTOL;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    void CreateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling) const
    {
        double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;

        if(     r_geom.GetGeometryFamily() == GeometryData::Kratos_Quadrilateral
            || (r_geom.GetGeometryFamily() == GeometryData::Kratos_NURBS && r_geom.GetGeometryType() == GeometryData::Kratos_Bezier2D) )
        {
            if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Quadrilateral)
            {
                Xmin = -1.0; Xmax = 1.0;
                Ymin = -1.0; Ymax = 1.0;
            }
            else
            {
                Xmin = 0.0; Xmax = 1.0;
                Ymin = 0.0; Ymax = 1.0;
            }

            double dX = (Xmax - Xmin) / nsampling;
            double dY = (Ymax - Ymin) / nsampling;

            SamplingPoints.reserve((nsampling+1) * (nsampling+1));
            CoordinatesArrayType X;
            PointType P;
            for(std::size_t i = 0; i < nsampling + 1; ++i)
            {
                X[0] = Xmin + i*dX;
                for(std::size_t j = 0; j < nsampling + 1; ++j)
                {
                    X[1] = Ymin + j*dY;
                    r_geom.GlobalCoordinates(P, X);
                    SamplingPoints.push_back(P);
                }
            }
        }
        else if( r_geom.GetGeometryFamily() == GeometryData::Kratos_Hexahedra
            ||  (r_geom.GetGeometryFamily() == GeometryData::Kratos_NURBS && r_geom.GetGeometryType() == GeometryData::Kratos_Bezier3D) )
        {
            if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Hexahedra)
            {
                Xmin = -1.0; Xmax = 1.0;
                Ymin = -1.0; Ymax = 1.0;
                Zmin = -1.0; Zmax = 1.0;
            }
            else
            {
                Xmin = 0.0; Xmax = 1.0;
                Ymin = 0.0; Ymax = 1.0;
                Zmin = 0.0; Zmax = 1.0;
            }

            double dX = (Xmax - Xmin) / nsampling;
            double dY = (Ymax - Ymin) / nsampling;
            double dZ = (Zmax - Zmin) / nsampling;

            SamplingPoints.reserve((nsampling+1) * (nsampling+1) * (nsampling+1));
            CoordinatesArrayType X;
            PointType P;
            for(std::size_t i = 0; i < nsampling + 1; ++i)
            {
                X[0] = Xmin + i*dX;
                for(std::size_t j = 0; j < nsampling + 1; ++j)
                {
                    X[1] = Ymin + j*dY;
                    for(std::size_t k = 0; k < nsampling + 1; ++k)
                    {
                        X[2] = Zmin + k*dZ;
                        r_geom.GlobalCoordinates(P, X);
                        SamplingPoints.push_back(P);
                    }
                }
            }
        }
    }


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
    BRep& operator=(BRep const& rOther);

    ///@}

}; // Class BRep

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BREP_H_INCLUDED  defined
