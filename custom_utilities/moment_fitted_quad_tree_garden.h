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
//  Date:            10 Mar 2017
//


#if !defined(KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_GARDEN_H_INCLUDED )
#define KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_GARDEN_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "custom_utilities/quadrature_utility.h"
#include "includes/model_part.h"
#include "includes/geometrical_object.h"
#include "utilities/math_utils.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_geometries/finite_cell_geometry.h"
#include "custom_utilities/quad_tree.h"


namespace Kratos
{

/// Short class definition.
/** A special type of quad tree garden. It uses quad-tree to compute the moment fitted quadrature of the element. For details see:
Hoang-Giang Bui, D. Schillinger, G. Meschke, Finite Cell Method for Plasticity using Moment-Fitted Quadrature Technique, in preparation.
*/
class MomentFittedQuadTreeGarden : public QuadTreeGarden
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  MomentFittedQuadTreeGarden
    KRATOS_CLASS_POINTER_DEFINITION(MomentFittedQuadTreeGarden);

    typedef QuadTreeGarden BaseType;

    typedef typename GeometricalObject::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor to construct the sub-cells based on Gauss quadrature.
    MomentFittedQuadTreeGarden(Element::Pointer& p_elem, const int& integration_order) : BaseType(p_elem)
    {
        if( p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral2D9
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
         || p_elem->GetGeometry().GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9 )
        {
            BaseType::ConstructSubCellsForQuadBasedOnGaussQuadrature(BaseType::mpTreeNodes, integration_order);
        }

        // check if the integration_points is inside the sub-cells
        GeometryData::IntegrationMethod ThisIntegrationMethod
                = LevelSet::GetIntegrationMethod(integration_order);

        const GeometryType::IntegrationPointsArrayType& integration_points
                = p_elem->GetGeometry().IntegrationPoints( ThisIntegrationMethod );

        for(std::size_t i = 0; i < BaseType::mpTreeNodes.size(); ++i)
        {
            if(!mpTreeNodes[i]->IsInside(integration_points[i]))
                KRATOS_THROW_ERROR(std::logic_error, "The integration_point is not inside the tree, error at", i)
        }
    }

    /// Destructor.
    virtual ~ MomentFittedQuadTreeGarden() {}


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MomentFittedQuadTreeGarden";
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

private:

    /// Assignment operator.
     MomentFittedQuadTreeGarden& operator=( MomentFittedQuadTreeGarden const& rOther);

    /// Copy constructor.
     MomentFittedQuadTreeGarden( MomentFittedQuadTreeGarden const& rOther);

}; // Class  MomentFittedQuadTreeGarden


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MomentFittedQuadTreeGarden& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDegree>
inline std::ostream& operator << (std::ostream& rOStream, const  MomentFittedQuadTreeGarden& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_FINITE_CELL_APPLICATION_MOMENT_FITTED_QUAD_TREE_GARDEN_H_INCLUDED defined
