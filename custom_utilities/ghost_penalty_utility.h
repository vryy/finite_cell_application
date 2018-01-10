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
//  Date:            6 Jan 2018
//


#if !defined(KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED )
#define  KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <limits>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>


// Project includes
#include "includes/define.h"
#include "includes/deprecated_variables.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "custom_algebra/brep.h"
#include "custom_conditions/ghost_penalty_condition.h"


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
/** Search and create the ghost penalty conditions
*/
class GhostPenaltyUtility
{
private:

    static const int msEdgesT3[][2];
    static const int msEdgesT6[][3];
    static const int msEdgesQ4[][2];
    static const int msEdgesQ8[][3];
    static const int msEdgesQ9[][3];
    static const int msFacesT4[][3];
    static const int msFacesT10[][6];
    static const int msFacesH8[][4];
    static const int msFacesH20[][8];
    static const int msFacesH27[][9];

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GhostPenaltyUtility
    KRATOS_CLASS_POINTER_DEFINITION(GhostPenaltyUtility);

    typedef Element::GeometryType GeometryType;

    typedef GeometryType::PointType NodeType;

    typedef NodeType::PointType PointType;

    typedef NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GhostPenaltyUtility() {}

    /// Destructor.
    virtual ~GhostPenaltyUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Setup the ghost penalty conditions for the model_part
    static ConditionsContainerType SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties);

    /// Setup the ghost penalty conditions for the model_part on a subset of elements
    static ConditionsContainerType SetUpSurfacePenaltyConditions(ModelPart& r_model_part,
        ModelPart::ElementsContainerType& pElements,
        GhostPenaltyCondition::Pointer p_sample_condition, const BRep& r_brep,
        std::size_t& lastCondId, Properties::Pointer pProperties);

    /// Setup the ghost penalty conditions between an element and its neighbour
    /// Note: user must call FindElementalNeighboursProcess(model.model_part, 2, 10).Execute() to setup first the neighbour elements
    static ConditionsContainerType SetUpSurfacePenaltyConditions(Element::Pointer p_element,
        GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties);

    /// Setup the ghost penalty conditions between two elements
    static Condition::Pointer SetUpSurfacePenaltyCondition(Element::Pointer p_element_1,
        Element::Pointer p_element_2, GhostPenaltyCondition::Pointer p_sample_condition,
        const BRep& r_brep, std::size_t& lastCondId, Properties::Pointer pProperties);

    /// Compute the shape function gradient in the normal direction, the edge geometry must be on the edge of the element
    static void ComputeShapeFunctionNormalGradient(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    /// Probe the neighbour elements of an element
    static void ProbeNeighbourElements(Element::Pointer p_element);


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
        return "Ghost Penalty Utility";
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


    /// Find the common edge between the two geometries, the index of the nodes on the edge are the local index w.r.t geometry 1
    static std::pair<GeometryData::KratosGeometryType, std::vector<std::size_t> > FindCommonEdge(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonEdgeT3(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonEdgeT6(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonEdgeQ4(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonEdgeQ8(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonFaceT4(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonFaceT10(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonFaceH8(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonFaceH20(Element& r_elem_1, Element& r_elem_2);

    static std::vector<std::size_t> FindCommonFaceH27(Element& r_elem_1, Element& r_elem_2);

    static bool IsBelongedToGeometry(const std::size_t& node_id, GeometryType& r_geom);

    static bool BuildMapEdgeNodeIndexToElementNodeIndex(std::map<std::size_t, std::size_t>& map_edge_node_index_to_element_node_index,
        GeometryType& r_edge_geometry, GeometryType& r_element_geometry, const int* nodes_on_edge);

    static void ComputeShapeFunctionNormalGradientT3(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientT6(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientQ4(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientQ8(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientQ9(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientT4(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientT10(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientH8(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientH20(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static void ComputeShapeFunctionNormalGradientH27(Matrix& dNdn, GeometryType& r_element_geometry,
        GeometryType& r_edge_geometry, const GeometryType::IntegrationPointsArrayType& integration_points);

    static bool IsSame(GeometryType& r_geom, const std::vector<std::size_t>& nodes);

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
    GhostPenaltyUtility& operator=(GhostPenaltyUtility const& rOther);

    /// Copy constructor.
    GhostPenaltyUtility(GhostPenaltyUtility const& rOther);


    ///@}

}; // Class GhostPenaltyUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, GhostPenaltyUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const GhostPenaltyUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GHOST_PENALTY_UTILITY_H_INCLUDED  defined
