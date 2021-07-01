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
//  Date:            16 Jun 2017
//


#if !defined(KRATOS_FINITE_CELL_GEOMETRY_UTILITY_H_INCLUDED )
#define  KRATOS_FINITE_CELL_GEOMETRY_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "custom_geometries/finite_cell_geometry.h"
#ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
#include "custom_geometries/finite_cell_geo_2d_bezier.h"
#include "custom_geometries/finite_cell_geo_2d_bezier_3.h"
#include "custom_geometries/finite_cell_geo_3d_bezier.h"
#endif

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
/** class for auxilliary routines
*/
class FiniteCellGeometryUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FiniteCellGeometryUtility
    KRATOS_CLASS_POINTER_DEFINITION(FiniteCellGeometryUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef typename GeometryType::IntegrationPointType IntegrationPointType;

    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef typename GeometryType::IntegrationPointsContainerType IntegrationPointsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FiniteCellGeometryUtility() {}

    /// Destructor.
    virtual ~FiniteCellGeometryUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Helper function to assign the geometry_data for finite_cell_geometry
    static void AssignGeometryData(GeometryType& r_geom,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const Vector& rWeights)
    {
        try
        {
            if(r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D3)
            {
                typedef FiniteCellGeometry<Triangle2D3<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D6)
            {
                typedef FiniteCellGeometry<Triangle2D6<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
            {
                typedef FiniteCellGeometry<Quadrilateral2D4<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8)
            {
                typedef FiniteCellGeometry<Quadrilateral2D8<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9)
            {
                typedef FiniteCellGeometry<Quadrilateral2D9<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
            {
                typedef FiniteCellGeometry<Tetrahedra3D4<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10)
            {
                typedef FiniteCellGeometry<Tetrahedra3D10<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
            {
                typedef FiniteCellGeometry<Hexahedra3D8<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D20)
            {
                typedef FiniteCellGeometry<Hexahedra3D20<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D27)
            {
                typedef FiniteCellGeometry<Hexahedra3D27<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, rWeights);
            }
            #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
            // here we don't support the Bezier geometry because they requires the full integration point definition -- hbui (it can be explained better but i can't recall for now'
            #endif
            else
            {
                std::stringstream ss;
                ss << "The geometry " << typeid(r_geom).name() << " of type " << r_geom.GetGeometryType() << " is not supported";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }
        catch(const std::bad_cast& e)
        {
            std::cout << "WARNING: the input geometry is not the FiniteCellGeometry. Hence the quadrature is not assigned" << std::endl;
        }
    }


    /// Helper function to assign the geometry_data for finite_cell_geometry
    static void AssignGeometryData(GeometryType& r_geom,
            const GeometryData::IntegrationMethod& ElementalIntegrationMethod,
            const IntegrationPointsArrayType& integration_points)
    {
        try
        {
            if(r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D3)
            {
                typedef FiniteCellGeometry<Triangle2D3<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Triangle2D6)
            {
                typedef FiniteCellGeometry<Triangle2D6<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
            {
                typedef FiniteCellGeometry<Quadrilateral2D4<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8)
            {
                typedef FiniteCellGeometry<Quadrilateral2D8<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9)
            {
                typedef FiniteCellGeometry<Quadrilateral2D9<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
            {
                typedef FiniteCellGeometry<Tetrahedra3D4<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D10)
            {
                typedef FiniteCellGeometry<Tetrahedra3D10<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
            {
                typedef FiniteCellGeometry<Hexahedra3D8<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D20)
            {
                typedef FiniteCellGeometry<Hexahedra3D20<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D27)
            {
                typedef FiniteCellGeometry<Hexahedra3D27<NodeType> > FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            #ifdef ENABLE_FINITE_CELL_ISOGEOMETRIC
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Bezier2D)
            {
                typedef FiniteCellGeo2dBezier<NodeType> FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Bezier2D3)
            {
                typedef FiniteCellGeo2dBezier3<NodeType> FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            else if(r_geom.GetGeometryType() == GeometryData::Kratos_Bezier3D)
            {
                typedef FiniteCellGeo3dBezier<NodeType> FiniteCellGeometryType;
                FiniteCellGeometryType& r_fc_geom = dynamic_cast<FiniteCellGeometryType&>(r_geom);
                r_fc_geom.AssignGeometryData(ElementalIntegrationMethod, integration_points);
            }
            #endif
            else
            {
                std::stringstream ss;
                ss << "The geometry " << typeid(r_geom).name() << " of type " << r_geom.GetGeometryType() << " is not supported";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }
        catch(const std::bad_cast& e)
        {
            std::cout << "WARNING: the input geometry is not the FiniteCellGeometry. Hence the quadrature is not assigned" << std::endl;
        }
    }

    /// Helper function to compute the global coordinates in the reference frame
    static CoordinatesArrayType& GlobalCoordinates0( const GeometryType& rGeometry, CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates )
    {
        if (rResult.size() != 3)
            rResult.resize(3, false);
        noalias( rResult ) = ZeroVector( 3 );

        Vector N( rGeometry.size() );
        rGeometry.ShapeFunctionsValues( N, LocalCoordinates );

        for ( std::size_t i = 0 ; i < rGeometry.size() ; ++i )
            noalias( rResult ) += N[i] * rGeometry[i].GetInitialPosition();

        return rResult;
    }

    /// Helper function to compute the global coordinates in the current frame
    static CoordinatesArrayType& GlobalCoordinates( const GeometryType& rGeometry, CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates )
    {
        if (rResult.size() != 3)
            rResult.resize(3, false);
        noalias( rResult ) = ZeroVector( 3 );

        Vector N( rGeometry.size() );
        rGeometry.ShapeFunctionsValues( N, LocalCoordinates );

        for ( std::size_t i = 0 ; i < rGeometry.size() ; ++i )
            noalias( rResult ) += N[i] * rGeometry[i];

        return rResult;
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
        return "Finite Cell Geometry Utility";
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
    FiniteCellGeometryUtility& operator=(FiniteCellGeometryUtility const& rOther);

    /// Copy constructor.
    FiniteCellGeometryUtility(FiniteCellGeometryUtility const& rOther);


    ///@}

}; // Class FiniteCellGeometryUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, FiniteCellGeometryUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FiniteCellGeometryUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_FINITE_CELL_GEOMETRY_UTILITY_H_INCLUDED  defined
